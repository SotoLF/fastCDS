#include "output_writer.hpp"
#include "utils.hpp"
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <ctime>
#include <iomanip>
#include <vector>

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace {

// Stock ofstreams use a ~8 KB stdio buffer, which means a million-row TSV
// hits write() ~hundreds of thousands of times. Bumping to 1 MiB cuts the
// syscall count by ~100×; on NVMe SSDs that's a noticeable wall-time win
// because each write() also entails a memcpy through the libc buffer.
//
// `pubsetbuf` must be called before the file is opened — that's why we
// can't just do `std::ofstream f(path)` and tweak after.
constexpr std::streamsize kWriteBufferBytes = 1 << 20;  // 1 MiB

class BufferedOfstream : public std::ofstream {
public:
    explicit BufferedOfstream(const std::string& path)
        : buf_(kWriteBufferBytes) {
        rdbuf()->pubsetbuf(buf_.data(), static_cast<std::streamsize>(buf_.size()));
        open(path);
    }

    // Without this, ~ofstream() runs AFTER buf_ destructs (derived members
    // first, then base) — by then the filebuf is pointing at freed memory
    // and the implicit flush writes garbage / nothing. Force the flush
    // here while buf_ is still alive.
    ~BufferedOfstream() override {
        if (is_open()) close();
    }
private:
    std::vector<char> buf_;
};

bool ensure_dir(const std::string& dir) {
    struct stat st;
    if (stat(dir.c_str(), &st) == 0) return S_ISDIR(st.st_mode);
    return mkdir(dir.c_str(), 0755) == 0;
}

std::string join(const std::string& dir, const std::string& name) {
    if (dir.empty()) return name;
    if (dir.back() == '/') return dir + name;
    return dir + "/" + name;
}

void w_int_or_na(std::ostream& f, uint32_t v) {
    if (v == 0) f << "NA"; else f << v;
}

void w_str_or_na(std::ostream& f, const std::string& s) {
    if (s.empty()) f << "NA"; else f << s;
}

void w_double(std::ostream& f, double v) {
    f << std::fixed << std::setprecision(4) << v;
}

std::string iso8601_now() {
    std::time_t t = std::time(nullptr);
    std::tm tm{};
#if defined(_WIN32)
    gmtime_s(&tm, &t);
#else
    gmtime_r(&t, &tm);
#endif
    char buf[32];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%SZ", &tm);
    return buf;
}

std::string output_kind_to_string(OutputKind k) {
    switch (k) {
        case OutputKind::CODING:  return "coding";
        case OutputKind::INTRONS: return "introns";
        case OutputKind::SPAN:    return "span";
        case OutputKind::ISOFORM: return "isoform";
        case OutputKind::BED12:   return "bed12";
        case OutputKind::ALL:     return "all";
    }
    return "unknown";
}

std::string bed_name(const DomainResult& r) {
    std::ostringstream oss;
    // Prefer the resolved protein id (set in summary) so ENST queries still
    // get a meaningful name; fall back to the raw input.
    const std::string& pid = r.summary.protein_id.empty()
        ? r.domain.protein_id : r.summary.protein_id;
    oss << pid;
    if (!r.domain.domain_id.empty()) oss << "_" << r.domain.domain_id;
    if (r.domain.has_domain()) oss << "_" << r.domain.start << "-" << r.domain.end;
    else oss << "_no_domain";
    return oss.str();
}

// ----- Headers ----------------------------------------------------------- //

const char* kSummaryHeader =
"input_id\tprotein_id\ttranscript_id\tgene_id\tgene_name\tdomain_id\tchrom\tstrand\t"
"aa_start\taa_end\tdomain_length_aa\tdomain_length_nt\tprotein_length_aa\t"
"domain_genomic_start\tdomain_genomic_end\tn_coding_segments\tfully_mapped\tno_domain_mode\t"
"input_id_type\tis_mane_select\tis_ensembl_canonical\t"
"cds_length_mismatch\tcds_nt_remainder\t"
"n_coding_exons_touched\tn_introns_spanned\tis_single_exon_domain\t"
"fraction_domain_in_largest_exon\tintron_burden_nt\tstatus\n";

// Header + columns shared by isoform_structure.tsv, domain_cds_segments.tsv, domain_introns.tsv.
const char* kFeatureTsvHeader =
"input_id\tgene_id\tgene_name\ttranscript_id\tprotein_id\tdomain_id\t"
"is_mane_select\tis_ensembl_canonical\tcds_length_mismatch\tcds_nt_remainder\t"
"chrom\tstrand\tfeature_type\tfeature_id\tfeature_part\texon_number\t"
"feature_genomic_start\tfeature_genomic_end\tfeature_length_nt\t"
"feature_order_genomic\tfeature_order_transcript\t"
"cds_nt_start\tcds_nt_end\taa_start_encoded\taa_end_encoded\t"
"overlaps_domain\t"
"domain_overlap_genomic_start\tdomain_overlap_genomic_end\t"
"domain_overlap_cds_nt_start\tdomain_overlap_cds_nt_end\t"
"domain_overlap_aa_start\tdomain_overlap_aa_end\t"
"domain_overlap_fraction_of_feature\tdomain_overlap_fraction_of_domain\t"
"plot_group\n";

const char* kUnmappedHeader =
"input_id\tprotein_id\taa_start\taa_end\tdomain_id\treason\n";

// ----- Per-result append helpers ----------------------------------------- //
// Each helper writes the rows produced by ONE DomainResult to the given
// stream. They are shared by both the one-shot write_all() path and the
// streaming StreamingWriter so that the on-disk format stays identical.

void append_summary_row(std::ostream& f, const DomainResult& r) {
    if (!r.mapped) {
        f << r.domain.input_id << "\t";
        w_str_or_na(f, r.domain.protein_id);
        f << "\tNA\tNA\tNA\t";
        w_str_or_na(f, r.domain.domain_id);
        f << "\tNA\tNA\t";
        if (r.domain.has_domain()) {
            f << r.domain.start << "\t" << r.domain.end
              << "\t" << (r.domain.end - r.domain.start + 1)
              << "\t" << ((r.domain.end - r.domain.start + 1) * 3);
        } else {
            f << "NA\tNA\tNA\tNA";
        }
        f << "\tNA\tNA\tNA\t0\tfalse\t"
          << (r.no_domain_mode ? "true" : "false")
          << "\tNA\tNA\tNA\tNA\tNA"
          // Phase 7 derived columns — NA for unmapped queries.
          << "\tNA\tNA\tNA\tNA\tNA"
          << "\t" << r.unmapped.reason << "\n";
        return;
    }
    const auto& s = r.summary;
    f << s.input_id << "\t";
    w_str_or_na(f, s.protein_id);
    f << "\t" << s.transcript_id << "\t";
    w_str_or_na(f, s.gene_id);
    f << "\t"; w_str_or_na(f, s.gene_name);
    f << "\t"; w_str_or_na(f, s.domain_id);
    f << "\t" << s.chrom << "\t" << s.strand << "\t";
    if (s.no_domain_mode) {
        f << "NA\tNA\tNA\tNA";
    } else {
        f << s.aa_start << "\t" << s.aa_end
          << "\t" << s.domain_length_aa << "\t" << s.domain_length_nt;
    }
    f << "\t" << s.protein_length_aa
      << "\t";
    if (s.no_domain_mode) {
        f << "NA\tNA\tNA";
    } else {
        f << s.domain_genomic_start << "\t" << s.domain_genomic_end
          << "\t" << s.n_coding_segments;
    }
    f << "\t" << (s.fully_mapped ? "true" : "false")
      << "\t" << (s.no_domain_mode ? "true" : "false")
      << "\t"; w_str_or_na(f, s.input_id_type);
    f << "\t" << tribool_to_string(s.is_mane_select)
      << "\t" << tribool_to_string(s.is_ensembl_canonical)
      << "\t" << (s.cds_length_mismatch ? "true" : "false")
      << "\t" << static_cast<int>(s.cds_nt_remainder)
      << "\t";
    // Phase 7 derived metrics. NA for no-domain queries (nothing to count).
    if (s.no_domain_mode) {
        f << "NA\tNA\tNA\tNA\tNA";
    } else {
        f << s.n_coding_exons_touched
          << "\t" << s.n_introns_spanned
          << "\t" << (s.is_single_exon_domain ? "true" : "false")
          << "\t"; w_double(f, s.fraction_domain_in_largest_exon);
        f << "\t" << s.intron_burden_nt;
    }
    f << "\t";
    const char* base = s.no_domain_mode ? "structure_only"
                                         : (s.fully_mapped ? "ok" : "partial");
    if (s.cds_length_mismatch && !s.no_domain_mode) {
        f << base << "_cds_mismatch";
    } else {
        f << base;
    }
    f << "\n";
}

void write_feature_tsv_row(std::ostream& f, const IsoformSegmentRow& s) {
    f << s.input_id << "\t";
    w_str_or_na(f, s.gene_id); f << "\t";
    w_str_or_na(f, s.gene_name); f << "\t";
    w_str_or_na(f, s.transcript_id); f << "\t";
    w_str_or_na(f, s.protein_id); f << "\t";
    w_str_or_na(f, s.domain_id); f << "\t"
      << tribool_to_string(s.is_mane_select) << "\t"
      << tribool_to_string(s.is_ensembl_canonical) << "\t"
      << (s.cds_length_mismatch ? "true" : "false") << "\t"
      << static_cast<int>(s.cds_nt_remainder) << "\t";
    f << s.chrom << "\t" << s.strand << "\t"
      << plot_feature_to_string(s.feature_type) << "\t"
      << s.feature_id << "\t"
      << s.feature_part << "\t";
    w_int_or_na(f, s.exon_number); f << "\t"
      << s.feature_genomic_start << "\t" << s.feature_genomic_end << "\t"
      << s.feature_length_nt << "\t"
      << s.feature_order_genomic << "\t" << s.feature_order_transcript << "\t";
    if (s.has_cds_coords) {
        f << s.cds_nt_start << "\t" << s.cds_nt_end << "\t"
          << s.aa_start_encoded << "\t" << s.aa_end_encoded;
    } else {
        f << "NA\tNA\tNA\tNA";
    }
    f << "\t";
    if (s.no_domain_mode) f << "NA";
    else f << overlap_kind_to_string(s.overlap);
    f << "\t";
    if (!s.no_domain_mode && s.has_overlap_coords) {
        f << s.domain_overlap_genomic_start << "\t" << s.domain_overlap_genomic_end << "\t"
          << s.domain_overlap_cds_nt_start << "\t" << s.domain_overlap_cds_nt_end << "\t"
          << s.domain_overlap_aa_start << "\t" << s.domain_overlap_aa_end << "\t";
        w_double(f, s.domain_overlap_fraction_of_feature); f << "\t";
        w_double(f, s.domain_overlap_fraction_of_domain);
    } else {
        f << "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA";
    }
    f << "\t" << s.plot_group << "\n";
}

void append_isoform_rows(std::ostream& f, const DomainResult& r) {
    if (!r.mapped) return;
    for (const auto& s : r.isoform_segments) write_feature_tsv_row(f, s);
}

void append_coding_rows(std::ostream& f, const DomainResult& r) {
    if (!r.mapped) return;
    for (const auto& s : r.isoform_segments) {
        if (s.feature_type == PlotFeatureType::CDS) write_feature_tsv_row(f, s);
    }
}

void append_introns_rows(std::ostream& f, const DomainResult& r) {
    if (!r.mapped) return;
    for (const auto& s : r.isoform_segments) {
        if (s.feature_type == PlotFeatureType::INTRON) write_feature_tsv_row(f, s);
    }
}

void append_coding_bed_rows(std::ostream& f, const DomainResult& r) {
    if (!r.mapped || r.no_domain_mode) return;
    std::string name = bed_name(r);
    for (const auto& s : r.isoform_segments) {
        if (s.feature_type != PlotFeatureType::CDS) continue;
        if (s.overlap != DomainOverlapKind::CODING_OVERLAP) continue;
        f << s.chrom << "\t" << (s.feature_genomic_start - 1) << "\t" << s.feature_genomic_end
          << "\t" << name << "\t0\t" << s.strand << "\n";
    }
}

void append_introns_bed_rows(std::ostream& f, const DomainResult& r) {
    if (!r.mapped || r.no_domain_mode) return;
    std::string name = bed_name(r);
    for (const auto& s : r.isoform_segments) {
        if (s.feature_type != PlotFeatureType::INTRON) continue;
        if (s.overlap != DomainOverlapKind::INSIDE_DOMAIN_GENOMIC_SPAN) continue;
        f << s.chrom << "\t" << (s.feature_genomic_start - 1) << "\t" << s.feature_genomic_end
          << "\t" << name << "\t0\t" << s.strand << "\n";
    }
}

void append_span_bed_row(std::ostream& f, const DomainResult& r) {
    if (!r.mapped || !r.has_span) return;
    std::string name = bed_name(r);
    const auto& s = r.span;
    f << s.chrom << "\t" << (s.genomic_start - 1) << "\t" << s.genomic_end
      << "\t" << name << "\t0\t" << s.strand << "\n";
}

void append_bed12_row(std::ostream& f, const DomainResult& r) {
    if (!r.mapped || r.no_domain_mode || !r.has_span) return;
    std::vector<std::pair<uint32_t, uint32_t>> blocks;
    for (const auto& s : r.isoform_segments) {
        if (s.feature_type != PlotFeatureType::CDS) continue;
        if (s.overlap != DomainOverlapKind::CODING_OVERLAP) continue;
        if (!s.has_overlap_coords) continue;
        blocks.emplace_back(s.domain_overlap_genomic_start,
                            s.domain_overlap_genomic_end);
    }
    if (blocks.empty()) return;
    std::sort(blocks.begin(), blocks.end());

    uint32_t chrom_start = r.span.genomic_start - 1; // 0-based
    uint32_t chrom_end   = r.span.genomic_end;       // half-open
    std::string name = bed_name(r);
    f << r.span.chrom << "\t" << chrom_start << "\t" << chrom_end << "\t"
      << name << "\t0\t" << r.span.strand << "\t"
      << chrom_start << "\t" << chrom_end << "\t"
      << "255,0,0\t" << blocks.size() << "\t";
    for (size_t i = 0; i < blocks.size(); ++i) {
        if (i) f << ",";
        f << (blocks[i].second - blocks[i].first + 1);
    }
    f << ",\t";
    for (size_t i = 0; i < blocks.size(); ++i) {
        if (i) f << ",";
        f << (blocks[i].first - 1 - chrom_start);
    }
    f << ",\n";
}

void append_unmapped_row(std::ostream& f, const DomainResult& r) {
    if (r.mapped) return;
    const auto& u = r.unmapped;
    f << u.input_id << "\t" << u.protein_id << "\t";
    if (r.domain.has_domain()) f << u.aa_start << "\t" << u.aa_end;
    else f << "NA\tNA";
    f << "\t";
    w_str_or_na(f, u.domain_id);
    f << "\t" << u.reason << "\n";
}

void write_metadata_file(std::ostream& f,
                         OutputKind kind,
                         size_t total, size_t mapped,
                         size_t unmapped, size_t no_domain,
                         const std::string& source,
                         const std::vector<std::string>& cli_args) {
    f << "{\n";
    f << "  \"tool\": \"prot2exon\",\n";
    f << "  \"version\": \"2.2.0\",\n";
    f << "  \"timestamp_utc\": \"" << iso8601_now() << "\",\n";
    f << "  \"output_kind\": \"" << output_kind_to_string(kind) << "\",\n";
    f << "  \"annotation_source\": \"" << source << "\",\n";
    f << "  \"index_format_version\": " << INDEX_FORMAT_VERSION << ",\n";
    f << "  \"coordinate_conventions\": {\n";
    f << "    \"bed\": \"0-based half-open\",\n";
    f << "    \"tsv\": \"1-based inclusive (genomic and CDS nt)\",\n";
    f << "    \"aa\": \"1-based inclusive\"\n";
    f << "  },\n";
    f << "  \"query_counts\": { \"total\": " << total
      << ", \"mapped\": " << mapped
      << ", \"unmapped\": " << unmapped
      << ", \"no_domain_mode\": " << no_domain << " },\n";
    f << "  \"cli\": [";
    for (size_t i = 0; i < cli_args.size(); ++i) {
        if (i) f << ", ";
        f << "\"";
        for (char c : cli_args[i]) { if (c == '\\' || c == '"') f << '\\'; f << c; }
        f << "\"";
    }
    f << "]\n";
    f << "}\n";
}

// ----- One-shot writers (delegate to the per-result helpers) ------------- //
// Each writer is a simple sequential loop. The caller's `parallel sections`
// in write_all() runs the 7+ files concurrently on different threads, which
// is where the speedup actually comes from at scale.
//
// (An earlier per-thread-buffer experiment was reverted: when called from
// inside `parallel sections`, OpenMP nested parallelism is off by default,
// so the inner parallel collapsed to single-threaded and added overhead
// without gain. See git history for the failed prototype.)

ErrorCode write_summary(const std::string& path, const std::vector<DomainResult>& results) {
    BufferedOfstream f(path);
    if (!f.is_open()) return ErrorCode::FILE_NOT_FOUND;
    f << kSummaryHeader;
    for (const auto& r : results) append_summary_row(f, r);
    std::cerr << "Wrote " << path << std::endl;
    return ErrorCode::SUCCESS;
}

ErrorCode write_isoform_tsv(const std::string& path, const std::vector<DomainResult>& results) {
    BufferedOfstream f(path);
    if (!f.is_open()) return ErrorCode::FILE_NOT_FOUND;
    f << kFeatureTsvHeader;
    for (const auto& r : results) append_isoform_rows(f, r);
    std::cerr << "Wrote " << path << std::endl;
    return ErrorCode::SUCCESS;
}

ErrorCode write_coding_tsv(const std::string& path, const std::vector<DomainResult>& results) {
    BufferedOfstream f(path);
    if (!f.is_open()) return ErrorCode::FILE_NOT_FOUND;
    f << kFeatureTsvHeader;
    for (const auto& r : results) append_coding_rows(f, r);
    std::cerr << "Wrote " << path << std::endl;
    return ErrorCode::SUCCESS;
}

ErrorCode write_introns_tsv(const std::string& path, const std::vector<DomainResult>& results) {
    BufferedOfstream f(path);
    if (!f.is_open()) return ErrorCode::FILE_NOT_FOUND;
    f << kFeatureTsvHeader;
    for (const auto& r : results) append_introns_rows(f, r);
    std::cerr << "Wrote " << path << std::endl;
    return ErrorCode::SUCCESS;
}

ErrorCode write_coding_bed(const std::string& path, const std::vector<DomainResult>& results) {
    BufferedOfstream f(path);
    if (!f.is_open()) return ErrorCode::FILE_NOT_FOUND;
    for (const auto& r : results) append_coding_bed_rows(f, r);
    std::cerr << "Wrote " << path << std::endl;
    return ErrorCode::SUCCESS;
}

ErrorCode write_introns_bed(const std::string& path, const std::vector<DomainResult>& results) {
    BufferedOfstream f(path);
    if (!f.is_open()) return ErrorCode::FILE_NOT_FOUND;
    for (const auto& r : results) append_introns_bed_rows(f, r);
    std::cerr << "Wrote " << path << std::endl;
    return ErrorCode::SUCCESS;
}

ErrorCode write_bed12(const std::string& path, const std::vector<DomainResult>& results) {
    BufferedOfstream f(path);
    if (!f.is_open()) return ErrorCode::FILE_NOT_FOUND;
    for (const auto& r : results) append_bed12_row(f, r);
    std::cerr << "Wrote " << path << std::endl;
    return ErrorCode::SUCCESS;
}

ErrorCode write_span_bed(const std::string& path, const std::vector<DomainResult>& results) {
    BufferedOfstream f(path);
    if (!f.is_open()) return ErrorCode::FILE_NOT_FOUND;
    for (const auto& r : results) append_span_bed_row(f, r);
    std::cerr << "Wrote " << path << std::endl;
    return ErrorCode::SUCCESS;
}

ErrorCode write_unmapped(const std::string& path, const std::vector<DomainResult>& results) {
    bool any = false;
    for (const auto& r : results) if (!r.mapped) { any = true; break; }
    if (!any) return ErrorCode::SUCCESS;
    BufferedOfstream f(path);
    if (!f.is_open()) return ErrorCode::FILE_NOT_FOUND;
    f << kUnmappedHeader;
    for (const auto& r : results) append_unmapped_row(f, r);
    std::cerr << "Wrote " << path << std::endl;
    return ErrorCode::SUCCESS;
}

ErrorCode write_metadata(const std::string& path,
                         OutputKind kind,
                         const std::vector<DomainResult>& results,
                         const std::string& source,
                         const std::vector<std::string>& cli_args) {
    BufferedOfstream f(path);
    if (!f.is_open()) return ErrorCode::FILE_NOT_FOUND;
    size_t mapped = 0, unmapped = 0, no_domain = 0;
    for (const auto& r : results) {
        if (r.mapped) ++mapped; else ++unmapped;
        if (r.no_domain_mode) ++no_domain;
    }
    write_metadata_file(f, kind, results.size(), mapped, unmapped, no_domain,
                        source, cli_args);
    std::cerr << "Wrote " << path << std::endl;
    return ErrorCode::SUCCESS;
}

} // namespace

namespace output {

ErrorCode write_all(const std::string& out_dir,
                    OutputKind kind,
                    const std::vector<DomainResult>& results,
                    const std::string& gtf_or_index_path,
                    const std::vector<std::string>& cli_args) {
    if (!ensure_dir(out_dir)) {
        std::cerr << "Error: cannot create output directory: " << out_dir << std::endl;
        return ErrorCode::FILE_NOT_FOUND;
    }

    // Each output file is an independent pass over results_, so write them
    // concurrently when OpenMP is available. When OMP is off, sections
    // collapse to sequential execution.
    ErrorCode rc_summary    = ErrorCode::SUCCESS;
    ErrorCode rc_coding_tsv = ErrorCode::SUCCESS;
    ErrorCode rc_coding_bed = ErrorCode::SUCCESS;
    ErrorCode rc_introns_tsv = ErrorCode::SUCCESS;
    ErrorCode rc_introns_bed = ErrorCode::SUCCESS;
    ErrorCode rc_span_bed   = ErrorCode::SUCCESS;
    ErrorCode rc_isoform    = ErrorCode::SUCCESS;
    ErrorCode rc_bed12      = ErrorCode::SUCCESS;
    ErrorCode rc_unmapped   = ErrorCode::SUCCESS;
    ErrorCode rc_metadata   = ErrorCode::SUCCESS;

    const bool want_coding  = (kind == OutputKind::CODING  || kind == OutputKind::ALL);
    const bool want_introns = (kind == OutputKind::INTRONS || kind == OutputKind::ALL);
    const bool want_span    = (kind == OutputKind::SPAN    || kind == OutputKind::ALL);
    const bool want_isoform = (kind == OutputKind::ISOFORM || kind == OutputKind::ALL);
    const bool want_bed12   = (kind == OutputKind::BED12   || kind == OutputKind::ALL);
    const bool want_meta    = (kind == OutputKind::ALL);

#ifdef USE_OPENMP
    #pragma omp parallel sections
#endif
    {
#ifdef USE_OPENMP
        #pragma omp section
#endif
        rc_summary = write_summary(join(out_dir, "domain_mapping_summary.tsv"), results);

#ifdef USE_OPENMP
        #pragma omp section
#endif
        if (want_coding)
            rc_coding_tsv = write_coding_tsv(join(out_dir, "domain_cds_segments.tsv"), results);

#ifdef USE_OPENMP
        #pragma omp section
#endif
        if (want_coding)
            rc_coding_bed = write_coding_bed(join(out_dir, "domain_cds_segments.bed"), results);

#ifdef USE_OPENMP
        #pragma omp section
#endif
        if (want_introns)
            rc_introns_tsv = write_introns_tsv(join(out_dir, "domain_introns.tsv"), results);

#ifdef USE_OPENMP
        #pragma omp section
#endif
        if (want_introns)
            rc_introns_bed = write_introns_bed(join(out_dir, "domain_introns.bed"), results);

#ifdef USE_OPENMP
        #pragma omp section
#endif
        if (want_span)
            rc_span_bed = write_span_bed(join(out_dir, "domain_span_with_introns.bed"), results);

#ifdef USE_OPENMP
        #pragma omp section
#endif
        if (want_isoform)
            rc_isoform = write_isoform_tsv(join(out_dir, "isoform_structure.tsv"), results);

#ifdef USE_OPENMP
        #pragma omp section
#endif
        if (want_bed12)
            rc_bed12 = write_bed12(join(out_dir, "domain_blocks.bed12"), results);

#ifdef USE_OPENMP
        #pragma omp section
#endif
        rc_unmapped = write_unmapped(join(out_dir, "unmapped_domains.tsv"), results);

#ifdef USE_OPENMP
        #pragma omp section
#endif
        if (want_meta)
            rc_metadata = write_metadata(join(out_dir, "run_metadata.json"), kind,
                                         results, gtf_or_index_path, cli_args);
    }

    for (ErrorCode rc : {rc_summary, rc_coding_tsv, rc_coding_bed,
                         rc_introns_tsv, rc_introns_bed, rc_span_bed,
                         rc_isoform, rc_bed12, rc_unmapped, rc_metadata}) {
        if (rc != ErrorCode::SUCCESS) return rc;
    }
    return ErrorCode::SUCCESS;
}

// --------------------------- StreamingWriter ----------------------------- //

StreamingWriter::StreamingWriter(std::string out_dir, OutputKind kind)
    : out_dir_(std::move(out_dir)), kind_(kind) {}

StreamingWriter::~StreamingWriter() = default;

ErrorCode StreamingWriter::open() {
    if (!ensure_dir(out_dir_)) {
        std::cerr << "Error: cannot create output directory: " << out_dir_ << std::endl;
        return ErrorCode::FILE_NOT_FOUND;
    }

    auto open_with_header = [&](std::ofstream& s, const char* name, const char* header) {
        s.open(join(out_dir_, name));
        if (!s.is_open()) return false;
        if (header) s << header;
        return true;
    };

    if (!open_with_header(summary_, "domain_mapping_summary.tsv", kSummaryHeader))
        return ErrorCode::FILE_NOT_FOUND;

    if (kind_ == OutputKind::CODING || kind_ == OutputKind::ALL) {
        if (!open_with_header(coding_tsv_, "domain_cds_segments.tsv", kFeatureTsvHeader))
            return ErrorCode::FILE_NOT_FOUND;
        if (!open_with_header(coding_bed_, "domain_cds_segments.bed", nullptr))
            return ErrorCode::FILE_NOT_FOUND;
    }
    if (kind_ == OutputKind::INTRONS || kind_ == OutputKind::ALL) {
        if (!open_with_header(introns_tsv_, "domain_introns.tsv", kFeatureTsvHeader))
            return ErrorCode::FILE_NOT_FOUND;
        if (!open_with_header(introns_bed_, "domain_introns.bed", nullptr))
            return ErrorCode::FILE_NOT_FOUND;
    }
    if (kind_ == OutputKind::SPAN || kind_ == OutputKind::ALL) {
        if (!open_with_header(span_bed_, "domain_span_with_introns.bed", nullptr))
            return ErrorCode::FILE_NOT_FOUND;
    }
    if (kind_ == OutputKind::ISOFORM || kind_ == OutputKind::ALL) {
        if (!open_with_header(isoform_tsv_, "isoform_structure.tsv", kFeatureTsvHeader))
            return ErrorCode::FILE_NOT_FOUND;
    }
    if (kind_ == OutputKind::BED12 || kind_ == OutputKind::ALL) {
        if (!open_with_header(bed12_, "domain_blocks.bed12", nullptr))
            return ErrorCode::FILE_NOT_FOUND;
    }

    opened_ = true;
    return ErrorCode::SUCCESS;
}

void StreamingWriter::ensure_unmapped_stream() {
    if (unmapped_) return;
    unmapped_ = std::make_unique<std::ofstream>(join(out_dir_, "unmapped_domains.tsv"));
    if (unmapped_->is_open()) {
        *unmapped_ << kUnmappedHeader;
    }
}

void StreamingWriter::append(const std::vector<DomainResult>& chunk) {
    for (const auto& r : chunk) {
        ++total_;
        if (r.mapped) ++mapped_; else ++n_unmapped_;
        if (r.no_domain_mode) ++no_domain_;

        append_summary_row(summary_, r);

        if (kind_ == OutputKind::CODING || kind_ == OutputKind::ALL) {
            append_coding_rows(coding_tsv_, r);
            append_coding_bed_rows(coding_bed_, r);
        }
        if (kind_ == OutputKind::INTRONS || kind_ == OutputKind::ALL) {
            append_introns_rows(introns_tsv_, r);
            append_introns_bed_rows(introns_bed_, r);
        }
        if (kind_ == OutputKind::SPAN || kind_ == OutputKind::ALL) {
            append_span_bed_row(span_bed_, r);
        }
        if (kind_ == OutputKind::ISOFORM || kind_ == OutputKind::ALL) {
            append_isoform_rows(isoform_tsv_, r);
        }
        if (kind_ == OutputKind::BED12 || kind_ == OutputKind::ALL) {
            append_bed12_row(bed12_, r);
        }
        if (!r.mapped) {
            ensure_unmapped_stream();
            append_unmapped_row(*unmapped_, r);
        }
    }
}

ErrorCode StreamingWriter::finalize(const std::string& gtf_or_index_path,
                                    const std::vector<std::string>& cli_args) {
    auto close_and_log = [&](std::ofstream& s, const char* name) {
        if (s.is_open()) {
            s.close();
            std::cerr << "Wrote " << join(out_dir_, name) << std::endl;
        }
    };
    close_and_log(summary_,     "domain_mapping_summary.tsv");
    close_and_log(coding_tsv_,  "domain_cds_segments.tsv");
    close_and_log(coding_bed_,  "domain_cds_segments.bed");
    close_and_log(introns_tsv_, "domain_introns.tsv");
    close_and_log(introns_bed_, "domain_introns.bed");
    close_and_log(span_bed_,    "domain_span_with_introns.bed");
    close_and_log(isoform_tsv_, "isoform_structure.tsv");
    close_and_log(bed12_,       "domain_blocks.bed12");
    if (unmapped_) close_and_log(*unmapped_, "unmapped_domains.tsv");

    if (kind_ == OutputKind::ALL) {
        std::string mpath = join(out_dir_, "run_metadata.json");
        std::ofstream m(mpath);
        if (!m.is_open()) return ErrorCode::FILE_NOT_FOUND;
        write_metadata_file(m, kind_, total_, mapped_, n_unmapped_, no_domain_,
                            gtf_or_index_path, cli_args);
        std::cerr << "Wrote " << mpath << std::endl;
    }
    return ErrorCode::SUCCESS;
}

} // namespace output
