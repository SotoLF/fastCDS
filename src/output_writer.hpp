#ifndef OUTPUT_WRITER_HPP
#define OUTPUT_WRITER_HPP

#include "domain_mapper.hpp"

#include <fstream>
#include <memory>

namespace output {

// Writes the requested set of files into out_dir based on output_kind.
// Always writes domain_mapping_summary.tsv and unmapped_domains.tsv (the latter
// only if there are unmapped rows). For ALL, also writes run_metadata.json.
ErrorCode write_all(const std::string& out_dir,
                    OutputKind kind,
                    const std::vector<DomainResult>& results,
                    const std::string& gtf_or_index_path,
                    const std::vector<std::string>& cli_args,
                    bool also_bed12 = false);

// Streaming writer: opens the requested set of files once with headers, then
// appends per-chunk rows so the caller can free each chunk's results before
// processing the next one. Use this to bound peak RAM at O(batch_size).
//
// Lifecycle:
//   1) open()           - create out_dir, open streams, write headers.
//   2) append(chunk)... - one call per chunk; rows are written, then chunk
//                         can be destroyed. Aggregates query counts.
//   3) finalize(...)    - close streams; for OutputKind::ALL also writes
//                         run_metadata.json with the accumulated counts.
//
// The unmapped_domains.tsv file is opened lazily on the first unmapped row
// (matching the non-streaming path: no file if no unmapped queries).
class StreamingWriter {
public:
    StreamingWriter(std::string out_dir, OutputKind kind, bool also_bed12 = false);
    ~StreamingWriter();

    StreamingWriter(const StreamingWriter&) = delete;
    StreamingWriter& operator=(const StreamingWriter&) = delete;

    ErrorCode open();
    void append(const std::vector<DomainResult>& chunk);
    ErrorCode finalize(const std::string& gtf_or_index_path,
                       const std::vector<std::string>& cli_args);

    size_t total_count() const     { return total_; }
    size_t mapped_count() const    { return mapped_; }
    size_t unmapped_count() const  { return n_unmapped_; }
    size_t no_domain_count() const { return no_domain_; }

private:
    void ensure_unmapped_stream();

    std::string out_dir_;
    OutputKind kind_;
    bool also_bed12_ = false;

    std::ofstream summary_;
    std::ofstream isoform_tsv_;
    std::ofstream coding_tsv_;
    std::ofstream introns_tsv_;
    std::ofstream coding_bed_;
    std::ofstream introns_bed_;
    std::ofstream span_bed_;
    std::ofstream bed12_;
    std::unique_ptr<std::ofstream> unmapped_;

    size_t total_ = 0;
    size_t mapped_ = 0;
    size_t n_unmapped_ = 0;
    size_t no_domain_ = 0;
    bool opened_ = false;
};

} // namespace output

#endif // OUTPUT_WRITER_HPP
