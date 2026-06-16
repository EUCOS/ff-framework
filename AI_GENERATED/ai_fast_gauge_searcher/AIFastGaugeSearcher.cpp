#include <algorithm>
#include <cerrno>
#include <csignal>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <future>
#include <iostream>
#include <list>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <BasisAlphaBuilder.h>
#include <BasisVector.h>
#include <Chunk.h>
#include <ChunkConsistencyChecker.h>
#include <GSOCoefficientMatrix.h>
#include <GSOCoefficientMatrixGenerator.h>
#include <LEEFT.h>
#include <MIBVGenerator.h>
#include <ModelBuilder.h>
#include <OutputWriter.h>

namespace {

volatile std::sig_atomic_t g_stop_requested = 0;

void HandleSignal(int) {
  g_stop_requested = 1;
}

struct Options {
  int dimension = 4;
  std::vector<int> orders;
  std::string output;
  uint64_t shard_index = 0;
  uint64_t shard_count = 1;
  unsigned threads = 1;
  std::string checkpoint;
  bool resume = false;
};

struct Checkpoint {
  uint64_t next_candidate = 0;
  uint64_t total_bvs = 0;
  uint64_t shard_consistent_bvs = 0;
  uint64_t consistent_models = 0;
  uint64_t unique_models = 0;
};

struct Task {
  uint64_t candidate = 0;
  std::vector<BasisVector> basis_vectors;
};

struct ModelRecord {
  LEEFT leeft;
  std::string output_block;
};

struct TaskResult {
  uint64_t candidate = 0;
  std::vector<ModelRecord> records;
  uint64_t consistent_models = 0;
};

bool StartsWith(const std::string& value, const std::string& prefix) {
  return value.compare(0, prefix.size(), prefix) == 0;
}

uint64_t ParseUInt64(const std::string& text, const std::string& name) {
  char* end = NULL;
  errno = 0;
  const unsigned long long value = std::strtoull(text.c_str(), &end, 10);
  if (errno != 0 || end == text.c_str() || *end != '\0') {
    throw std::runtime_error("Invalid integer for " + name + ": " + text);
  }
  return static_cast<uint64_t>(value);
}

int ParseInt(const std::string& text, const std::string& name) {
  const uint64_t value = ParseUInt64(text, name);
  if (value > static_cast<uint64_t>(2147483647)) {
    throw std::runtime_error("Integer too large for " + name + ": " + text);
  }
  return static_cast<int>(value);
}

std::vector<int> ParseOrders(const std::string& text) {
  std::vector<int> orders;
  std::stringstream ss(text);
  std::string item;
  while (std::getline(ss, item, ',')) {
    if (!item.empty()) {
      orders.push_back(ParseInt(item, "--orders"));
    }
  }
  if (orders.empty()) {
    throw std::runtime_error("--orders must contain at least one order");
  }
  return orders;
}

Options ParseOptions(int argc, char* argv[]) {
  Options options;
  for (int i = 1; i < argc; ++i) {
    const std::string arg(argv[i]);
    if (arg == "--dimension" && i + 1 < argc) {
      options.dimension = ParseInt(argv[++i], "--dimension");
    } else if (arg == "--orders" && i + 1 < argc) {
      options.orders = ParseOrders(argv[++i]);
    } else if (arg == "--output" && i + 1 < argc) {
      options.output = argv[++i];
    } else if (arg == "--shard-index" && i + 1 < argc) {
      options.shard_index = ParseUInt64(argv[++i], "--shard-index");
    } else if (arg == "--shard-count" && i + 1 < argc) {
      options.shard_count = ParseUInt64(argv[++i], "--shard-count");
    } else if (arg == "--threads" && i + 1 < argc) {
      options.threads = static_cast<unsigned>(ParseUInt64(argv[++i], "--threads"));
    } else if (arg == "--checkpoint" && i + 1 < argc) {
      options.checkpoint = argv[++i];
    } else if (arg == "--resume") {
      options.resume = true;
    } else if (arg == "--help") {
      std::cout
          << "Usage: AIFastGaugeSearcher --dimension 4 --orders 30 "
          << "--output <path> --shard-index <i> --shard-count <n> "
          << "--threads <n> --checkpoint <path> [--resume]\n";
      std::exit(0);
    } else {
      throw std::runtime_error("Unknown or incomplete argument: " + arg);
    }
  }

  if (options.dimension != 4) {
    throw std::runtime_error("AIFastGaugeSearcher currently supports only --dimension 4");
  }
  if (options.orders.empty()) {
    throw std::runtime_error("--orders is required");
  }
  if (options.output.empty()) {
    throw std::runtime_error("--output is required");
  }
  if (options.shard_count == 0) {
    throw std::runtime_error("--shard-count must be greater than zero");
  }
  if (options.shard_index >= options.shard_count) {
    throw std::runtime_error("--shard-index must be less than --shard-count");
  }
  if (options.threads == 0) {
    throw std::runtime_error("--threads must be greater than zero");
  }
  if (options.checkpoint.empty()) {
    options.checkpoint = options.output + ".checkpoint";
  }
  return options;
}

ModelBuilder BuildInitialConditions(int dimension) {
  ModelBuilder initial(dimension);
  initial.Load_S_Vector(dimension);
  initial.Load_Default_k_ij();
  GSOCoefficientMatrix initial_kij = initial.FFHS_Model().k_ij();
  if (!initial.Check_Model_Consistency()) {
    throw std::runtime_error("Bad initial model conditions");
  }
  initial.rFFHS_Model().Set_k_ij(initial_kij);
  return initial;
}

void FormKijExtensions(
    const std::vector<std::list<std::vector<int> > >& layer_extensions,
    int layer,
    std::vector<std::vector<int> >& current,
    std::list<std::vector<std::vector<int> > >& output) {
  if (layer == static_cast<int>(layer_extensions.size())) {
    output.push_back(current);
    return;
  }
  for (std::list<std::vector<int> >::const_iterator it =
           layer_extensions.at(layer).begin();
       it != layer_extensions.at(layer).end(); ++it) {
    current.at(layer) = *it;
    FormKijExtensions(layer_extensions, layer + 1, current, output);
  }
}

std::list<std::vector<std::vector<int> > > BuildKijExtensions(
    const Model& initial_model,
    const std::vector<int>& extension_orders) {
  std::vector<int> orders;
  for (int i = 0; i < static_cast<int>(initial_model.BV_Set().size()); ++i) {
    orders.push_back(initial_model.BV_Set().at(i).Order());
  }

  std::vector<std::list<std::vector<int> > > layer_extensions;
  for (int layer = 0; layer < static_cast<int>(extension_orders.size());
       ++layer) {
    orders.push_back(extension_orders.at(layer));
    GSOCoefficientMatrixGenerator generator(orders, layer);
    generator.Build_GSO_Coefficient_Extensions();
    layer_extensions.push_back(generator.GSO_Coefficient_Matrix_Extensions());
  }

  std::list<std::vector<std::vector<int> > > result;
  std::vector<std::vector<int> > current(extension_orders.size());
  FormKijExtensions(layer_extensions, 0, current, result);
  return result;
}

void AppendChunk(std::vector<int>& bv, const Chunk& chunk) {
  const std::vector<int>& values = chunk.BV_Chunk();
  bv.insert(bv.end(), values.begin(), values.end());
}

std::string RenderModelBlock(const Model& model,
                             const std::vector<std::vector<int> >& kij_extension,
                             int extensions) {
  std::ostringstream out;
  for (int a = 0; a < static_cast<int>(model.Gauge_Groups().size()); ++a) {
    out << model.Gauge_Groups().at(a).Name().Class() << " "
        << model.Gauge_Groups().at(a).Name().Rank() << " "
        << model.Gauge_Groups().at(a).Name().KM_Level() << " ";
  }
  out << "\n";

  for (int a = 0; a < static_cast<int>(model.MatterRepresentations().size()); ++a) {
    out << model.MatterRepresentations().at(a).Duplicates() << ": ";
    for (int b = 0;
         b < static_cast<int>(
                 model.MatterRepresentations().at(a).Rep_Dimension().size());
         ++b) {
      out << model.MatterRepresentations().at(a).Rep_Dimension().at(b).Dimension()
          << model.MatterRepresentations().at(a).Rep_Dimension().at(b).Triality()
          << " ";
    }
    out << "\n";
  }

  out << "U(1)'s: " << model.U1_Factors() << "\n";
  out << "ST SUSY: " << model.SUSY_States().size() << "\n";
  out << "\n";

  const int finish = static_cast<int>(model.BV_Set().size());
  const int start = finish - extensions;
  const int bv_size = static_cast<int>(model.BV_Set().at(0).BV().size());
  for (int a = start; a < finish; ++a) {
    for (int b = 0; b < bv_size; ++b) {
      out << model.BV_Set().at(a).BV().at(b) << " ";
    }
    out << "\n";
  }

  for (int a = 0; a < static_cast<int>(kij_extension.size()); ++a) {
    for (int b = 0; b < static_cast<int>(kij_extension.at(a).size()); ++b) {
      out << kij_extension.at(a).at(b) << " ";
    }
    out << "\n";
  }
  out << "\n";
  return out.str();
}

bool FileExistsAndNonEmpty(const std::string& path) {
  std::ifstream in(path.c_str());
  return in.good() && in.peek() != std::ifstream::traits_type::eof();
}

Checkpoint ReadCheckpoint(const std::string& path) {
  std::ifstream in(path.c_str());
  if (!in) {
    throw std::runtime_error("Could not open checkpoint for resume: " + path);
  }

  std::string version;
  in >> version;
  if (version != "AI_FAST_GAUGE_CHECKPOINT_V1") {
    throw std::runtime_error("Unrecognized checkpoint format: " + path);
  }

  Checkpoint checkpoint;
  std::string key;
  while (in >> key) {
    uint64_t value = 0;
    in >> value;
    if (key == "next_candidate") {
      checkpoint.next_candidate = value;
    } else if (key == "total_bvs") {
      checkpoint.total_bvs = value;
    } else if (key == "shard_consistent_bvs") {
      checkpoint.shard_consistent_bvs = value;
    } else if (key == "consistent_models") {
      checkpoint.consistent_models = value;
    } else if (key == "unique_models") {
      checkpoint.unique_models = value;
    }
  }
  return checkpoint;
}

void WriteCheckpoint(const std::string& path,
                     const Options& options,
                     const Checkpoint& checkpoint) {
  const std::string tmp_path = path + ".tmp";
  std::ofstream out(tmp_path.c_str(), std::ios::out | std::ios::trunc);
  if (!out) {
    throw std::runtime_error("Could not write checkpoint: " + tmp_path);
  }
  out << "AI_FAST_GAUGE_CHECKPOINT_V1\n";
  out << "dimension " << options.dimension << "\n";
  out << "shard_index " << options.shard_index << "\n";
  out << "shard_count " << options.shard_count << "\n";
  out << "threads " << options.threads << "\n";
  out << "next_candidate " << checkpoint.next_candidate << "\n";
  out << "total_bvs " << checkpoint.total_bvs << "\n";
  out << "shard_consistent_bvs " << checkpoint.shard_consistent_bvs << "\n";
  out << "consistent_models " << checkpoint.consistent_models << "\n";
  out << "unique_models " << checkpoint.unique_models << "\n";
  out << "timestamp " << static_cast<uint64_t>(std::time(NULL)) << "\n";
  out.close();
  if (std::rename(tmp_path.c_str(), path.c_str()) != 0) {
    throw std::runtime_error("Could not atomically replace checkpoint: " + path);
  }
}

TaskResult EvaluateTask(
    const Task& task,
    const Model& initial_model,
    int dimension,
    const std::list<std::vector<std::vector<int> > >& kij_extensions) {
  TaskResult result;
  result.candidate = task.candidate;

  bool linear_independence_checked = false;
  for (std::list<std::vector<std::vector<int> > >::const_iterator it =
           kij_extensions.begin();
       it != kij_extensions.end(); ++it) {
    ModelBuilder builder(dimension);
    for (int b = 1; b < static_cast<int>(initial_model.BV_Set().size()); ++b) {
      builder.Load_Basis_Vector(initial_model.BV_Set().at(b));
    }
    for (int b = 0; b < static_cast<int>(task.basis_vectors.size()); ++b) {
      builder.Load_Basis_Vector(task.basis_vectors.at(b));
    }

    if (!linear_independence_checked) {
      if (!builder.Check_Linear_Independence()) {
        break;
      }
      linear_independence_checked = true;
    }

    for (int b = 0;
         b < static_cast<int>(initial_model.k_ij().Numerators().size()); ++b) {
      builder.Load_k_ij_Row(initial_model.k_ij().Numerators().at(b));
    }
    for (int b = 0; b < static_cast<int>(it->size()); ++b) {
      builder.Load_k_ij_Row(it->at(b));
    }

    if (builder.Check_k_ij_Consistency()) {
      // Exactness note: all physics-sensitive work remains in the framework:
      // complete GSO/GGSO consistency, gauge group construction and naming, and
      // LEEFT creation. This driver only changes enumeration, batching, and I/O.
      builder.Build_Gauge_Group_Model();
      ModelRecord record;
      record.leeft = LEEFT(builder.FFHS_Model());
      record.output_block =
          RenderModelBlock(builder.FFHS_Model(), *it, task.basis_vectors.size());
      result.records.push_back(record);
      ++result.consistent_models;
    }
  }
  return result;
}

class AIFastGaugeScan {
 public:
  AIFastGaugeScan(const Options& options, const ModelBuilder& initial)
      : options_(options),
        initial_model_(initial.FFHS_Model()),
        initial_common_basis_alphas_(initial.Common_Basis_Alphas()),
        kij_extensions_(BuildKijExtensions(initial.FFHS_Model(), options.orders)),
        basis_vectors_(options.orders.size()) {
    if (options_.resume) {
      checkpoint_ = ReadCheckpoint(options_.checkpoint);
    }
  }

  int Run() {
    const bool append = options_.resume && FileExistsAndNonEmpty(options_.output);
    if (append) {
      RebuildUniqueLeeftsFromOutput();
    }

    char model_buffer[1 << 20];
    char md_buffer[1 << 16];
    model_out_.rdbuf()->pubsetbuf(model_buffer, sizeof(model_buffer));
    md_out_.rdbuf()->pubsetbuf(md_buffer, sizeof(md_buffer));
    model_out_.open(options_.output.c_str(),
                    append ? (std::ios::out | std::ios::app)
                           : (std::ios::out | std::ios::trunc));
    md_out_.open((options_.output + "MD").c_str(),
                 append ? (std::ios::out | std::ios::app)
                        : (std::ios::out | std::ios::trunc));
    if (!model_out_ || !md_out_) {
      throw std::runtime_error("Could not open output path: " + options_.output);
    }

    if (!append) {
      writer_.Write_BV_Search_Front_Info(
          model_out_, initial_model_, options_.dimension, options_.orders);
    }

    const std::clock_t start = std::clock();
    BuildExtensions(0);
    FlushBatch(global_candidate_);
    const double total_time =
        (std::clock() - start) / static_cast<double>(CLOCKS_PER_SEC);

    if (g_stop_requested) {
      FlushAndCheckpoint();
      std::cout << "Stop requested; checkpoint written at candidate "
                << checkpoint_.next_candidate << std::endl;
      return 2;
    }

    WriteEndMetadata(total_time);
    FlushAndCheckpoint();
    return 0;
  }

 private:
  Options options_;
  Model initial_model_;
  std::vector<BasisAlpha> initial_common_basis_alphas_;
  std::list<std::vector<std::vector<int> > > kij_extensions_;
  std::vector<BasisVector> basis_vectors_;
  OutputWriter writer_;
  std::ofstream model_out_;
  std::ofstream md_out_;
  std::set<LEEFT> unique_leefts_;
  Checkpoint checkpoint_;
  std::vector<Task> pending_tasks_;
  uint64_t global_candidate_ = 0;
  std::time_t last_checkpoint_time_ = 0;

  void RebuildUniqueLeeftsFromOutput() {
    std::ifstream in(options_.output.c_str());
    if (!in) {
      throw std::runtime_error("Could not read existing output for resume: " +
                               options_.output);
    }

    std::string line;
    bool in_models = false;
    while (std::getline(in, line)) {
      if (!in_models) {
        if (StartsWith(line, "--")) {
          in_models = true;
        }
        continue;
      }
      if (line.empty()) {
        continue;
      }
      if (StartsWith(line, "Total unique models:")) {
        break;
      }

      std::vector<std::string> leeft_data;
      leeft_data.push_back(line);
      while (std::getline(in, line)) {
        leeft_data.push_back(line);
        if (StartsWith(line, "ST SUSY:")) {
          break;
        }
      }
      if (leeft_data.size() >= 3 && StartsWith(leeft_data.back(), "ST SUSY:")) {
        unique_leefts_.insert(LEEFT(leeft_data));
      }

      const int skip_lines = 2 * static_cast<int>(options_.orders.size()) + 2;
      for (int i = 0; i < skip_lines && std::getline(in, line); ++i) {
      }
    }
  }

  void BuildExtensions(int layer) {
    if (g_stop_requested) {
      return;
    }
    if (layer < static_cast<int>(options_.orders.size())) {
      BuildBasisVectors(layer);
    } else {
      QueueCandidate();
    }
  }

  void BuildBasisVectors(int layer) {
    MIBVGenerator generator(options_.orders.at(layer), options_.dimension);
    if (layer == 0) {
      generator.Build_Gauge_Chunks(initial_common_basis_alphas_);
    } else {
      std::vector<BasisVector> all_basis_vectors = initial_model_.BV_Set();
      all_basis_vectors.reserve(all_basis_vectors.size() + layer);
      for (int i = 0; i < layer; ++i) {
        all_basis_vectors.push_back(basis_vectors_.at(i));
      }
      BasisAlphaBuilder ba_builder(all_basis_vectors);
      ba_builder.Build_Basis_Alphas();
      ba_builder.Build_Common_Basis_Alphas();
      generator.Build_Gauge_Chunks(ba_builder.Common_Basis_Alphas());
    }

    const std::list<Chunk>& lms = generator.SP_LMs();
    const std::list<Chunk>& comp = generator.SP_RMs_Compact();
    const std::list<Chunk>& obs = generator.RMs_Observable();
    const std::list<Chunk>& hid = generator.RMs_Hidden();
    if (lms.empty()) {
      return;
    }

    if (checkpoint_.total_bvs == 0 && layer == 0) {
      checkpoint_.total_bvs = static_cast<uint64_t>(lms.size()) *
                              static_cast<uint64_t>(obs.size()) *
                              static_cast<uint64_t>(comp.size()) *
                              static_cast<uint64_t>(hid.size());
    }

    ChunkConsistencyChecker checker(
        options_.orders.at(layer),
        initial_common_basis_alphas_.at(0).Denominator());

    const Chunk& lm = lms.front();
    for (std::list<Chunk>::const_iterator it_comp = comp.begin();
         it_comp != comp.end() && !g_stop_requested; ++it_comp) {
      for (std::list<Chunk>::const_iterator it_obs = obs.begin();
           it_obs != obs.end() && !g_stop_requested; ++it_obs) {
        for (std::list<Chunk>::const_iterator it_hid = hid.begin();
             it_hid != hid.end() && !g_stop_requested; ++it_hid) {
          if (!checker.Check_Modular_Invariance(lm, *it_obs, *it_comp, *it_hid)) {
            continue;
          }

          std::vector<int> new_bv;
          new_bv.reserve(lm.BV_Chunk().size() + it_obs->BV_Chunk().size() +
                         it_comp->BV_Chunk().size() + it_hid->BV_Chunk().size());
          AppendChunk(new_bv, lm);
          AppendChunk(new_bv, *it_obs);
          AppendChunk(new_bv, *it_comp);
          AppendChunk(new_bv, *it_hid);

          basis_vectors_.at(layer) =
              BasisVector(new_bv, options_.orders.at(layer), options_.dimension);
          BuildExtensions(layer + 1);
        }
      }
    }
  }

  void QueueCandidate() {
    const uint64_t candidate = global_candidate_++;
    if (candidate < checkpoint_.next_candidate) {
      return;
    }

    if (candidate % options_.shard_count != options_.shard_index) {
      if (pending_tasks_.empty()) {
        checkpoint_.next_candidate = candidate + 1;
        MaybeCheckpoint();
      }
      return;
    }

    Task task;
    task.candidate = candidate;
    task.basis_vectors = basis_vectors_;
    pending_tasks_.push_back(task);

    const size_t batch_size = std::max<size_t>(1, options_.threads * 2);
    if (pending_tasks_.size() >= batch_size) {
      FlushBatch(candidate + 1);
    }
  }

  void FlushBatch(uint64_t completed_through_candidate) {
    if (pending_tasks_.empty()) {
      return;
    }

    std::vector<std::future<TaskResult> > futures;
    futures.reserve(pending_tasks_.size());
    for (size_t i = 0; i < pending_tasks_.size(); ++i) {
      futures.push_back(std::async(
          std::launch::async,
          EvaluateTask,
          pending_tasks_.at(i),
          initial_model_,
          options_.dimension,
          kij_extensions_));

      if (futures.size() == options_.threads || i + 1 == pending_tasks_.size()) {
        for (size_t j = 0; j < futures.size(); ++j) {
          ApplyResult(futures.at(j).get());
        }
        futures.clear();
      }
    }

    checkpoint_.next_candidate = completed_through_candidate;
    pending_tasks_.clear();
    MaybeCheckpoint();
  }

  void ApplyResult(const TaskResult& result) {
    ++checkpoint_.shard_consistent_bvs;
    checkpoint_.consistent_models += result.consistent_models;
    for (size_t i = 0; i < result.records.size(); ++i) {
      const ModelRecord& record = result.records.at(i);
      const size_t previous_unique_count = unique_leefts_.size();
      unique_leefts_.insert(record.leeft);
      if (unique_leefts_.size() > previous_unique_count) {
        checkpoint_.unique_models = unique_leefts_.size();
        md_out_ << checkpoint_.consistent_models << " "
                << unique_leefts_.size() << std::endl;
        model_out_ << record.output_block;
      }
    }
  }

  void MaybeCheckpoint() {
    const std::time_t now = std::time(NULL);
    if (last_checkpoint_time_ == 0 || now - last_checkpoint_time_ >= 60) {
      FlushAndCheckpoint();
      last_checkpoint_time_ = now;
    }
  }

  void FlushAndCheckpoint() {
    checkpoint_.unique_models = unique_leefts_.size();
    model_out_.flush();
    md_out_.flush();
    WriteCheckpoint(options_.checkpoint, options_, checkpoint_);
  }

  void WriteEndMetadata(double total_time) {
    model_out_ << "Total unique models: " << unique_leefts_.size() << std::endl;
    model_out_ << "Total consistent models: " << checkpoint_.consistent_models
               << std::endl;
    model_out_ << "Total consistent BVs: " << checkpoint_.shard_consistent_bvs
               << std::endl;
    model_out_ << "Total BVs tested: " << checkpoint_.total_bvs << std::endl;
    const double models_per_sec =
        total_time > 0.0
            ? static_cast<double>(checkpoint_.consistent_models) / total_time
            : 0.0;
    model_out_ << "Total time: " << total_time << std::endl;
    model_out_ << "Models built per sec: " << models_per_sec << std::endl;
  }
};

}  // namespace

int main(int argc, char* argv[]) {
  try {
    const Options options = ParseOptions(argc, argv);
    std::signal(SIGTERM, HandleSignal);
    std::signal(SIGINT, HandleSignal);

    AIFastGaugeScan scan(options, BuildInitialConditions(options.dimension));
    return scan.Run();
  } catch (const std::exception& ex) {
    std::cerr << "AIFastGaugeSearcher error: " << ex.what() << std::endl;
    return 1;
  }
}
