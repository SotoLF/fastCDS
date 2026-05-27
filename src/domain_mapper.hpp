#ifndef DOMAIN_MAPPER_HPP
#define DOMAIN_MAPPER_HPP

#include "common.hpp"
#include "gtf_parser.hpp"
#include <functional>

struct ProteinDomain {
    std::string input_id;
    std::string protein_id;
    uint32_t start;              // 1-based aa; 0 means "no domain"
    uint32_t end;                // 1-based aa, inclusive
    std::string domain_id;

    ProteinDomain() : start(0), end(0) {}
    bool has_domain() const { return start != 0 && end != 0; }
};

// One row per input domain capturing the full genomic envelope (incl. introns).
struct SpanRow {
    std::string chrom;
    uint32_t genomic_start;
    uint32_t genomic_end;
    char strand;
    std::string input_id;
    std::string protein_id;
    std::string domain_id;
    std::string gene_name;
};

// One row per structural feature of the transcript.
// CDS that overlap the domain partially are split into multiple rows that
// SHARE the same feature_id and differ by feature_part.
struct IsoformSegmentRow {
    std::string input_id;
    std::string gene_id;
    std::string gene_name;
    std::string transcript_id;
    std::string protein_id;
    std::string domain_id;

    // Transcript-level flags (same value on every row of the same query).
    TriBool is_mane_select = TriBool::NA;
    TriBool is_ensembl_canonical = TriBool::NA;
    bool cds_length_mismatch = false;
    uint8_t cds_nt_remainder = 0; // 0, 1 or 2

    std::string chrom;
    char strand;

    PlotFeatureType feature_type;
    std::string feature_id;       // e.g. CDS_3 — stable across splits of the same GTF CDS
    uint32_t feature_part;        // 1..K when CDS is split by overlap; else 1
    uint32_t source_cds_index;    // 0-based index in translation order (CDS only); 0 for UTR/intron
    uint32_t exon_number;         // 0 -> NA
    uint32_t feature_genomic_start;
    uint32_t feature_genomic_end;
    uint32_t feature_length_nt;
    uint32_t feature_order_genomic;
    uint32_t feature_order_transcript;

    bool has_cds_coords;
    uint32_t cds_nt_start;
    uint32_t cds_nt_end;
    uint32_t aa_start_encoded;
    uint32_t aa_end_encoded;

    DomainOverlapKind overlap;
    bool has_overlap_coords;
    uint32_t domain_overlap_genomic_start;
    uint32_t domain_overlap_genomic_end;
    uint32_t domain_overlap_cds_nt_start;
    uint32_t domain_overlap_cds_nt_end;
    uint32_t domain_overlap_aa_start;
    uint32_t domain_overlap_aa_end;
    double domain_overlap_fraction_of_feature;
    double domain_overlap_fraction_of_domain;

    bool no_domain_mode;          // true when the input row has no domain at all
    std::string plot_group;

    IsoformSegmentRow()
        : strand('+'), feature_type(PlotFeatureType::CDS),
          feature_part(1), source_cds_index(0), exon_number(0),
          feature_genomic_start(0), feature_genomic_end(0), feature_length_nt(0),
          feature_order_genomic(0), feature_order_transcript(0),
          has_cds_coords(false), cds_nt_start(0), cds_nt_end(0),
          aa_start_encoded(0), aa_end_encoded(0),
          overlap(DomainOverlapKind::NO), has_overlap_coords(false),
          domain_overlap_genomic_start(0), domain_overlap_genomic_end(0),
          domain_overlap_cds_nt_start(0), domain_overlap_cds_nt_end(0),
          domain_overlap_aa_start(0), domain_overlap_aa_end(0),
          domain_overlap_fraction_of_feature(0.0),
          domain_overlap_fraction_of_domain(0.0),
          no_domain_mode(false) {}
};

struct MappingSummaryRow {
    std::string input_id;
    std::string protein_id;       // "" => no protein for this transcript
    std::string transcript_id;
    std::string gene_id;
    std::string gene_name;
    std::string domain_id;
    std::string chrom;
    char strand;
    uint32_t aa_start;            // 0 => no_domain
    uint32_t aa_end;
    uint32_t domain_length_aa;
    uint32_t domain_length_nt;
    uint32_t protein_length_aa;
    uint32_t domain_genomic_start;
    uint32_t domain_genomic_end;
    uint32_t n_coding_segments;
    bool fully_mapped;
    bool no_domain_mode;

    // Derived architecture metrics (zero/false in no_domain or unmapped rows).
    // Cheap to compute alongside the mapping; useful for atlas-scale analyses
    // so users don't have to re-aggregate from the segment TSVs every time.
    uint32_t n_coding_exons_touched = 0;        // distinct GTF CDS exons that overlap the domain
    uint32_t n_introns_spanned      = 0;        // introns between two coding_overlap CDS rows
    bool     is_single_exon_domain  = false;    // n_coding_exons_touched == 1
    double   fraction_domain_in_largest_exon = 0.0;  // [0,1]; sums of overlap_fraction_of_domain
    uint32_t intron_burden_nt       = 0;        // sum of lengths of the introns spanned

    // New in v2.2: transcript-level flags & CDS-mismatch.
    TriBool is_mane_select = TriBool::NA;
    TriBool is_ensembl_canonical = TriBool::NA;
    bool cds_length_mismatch = false;
    uint8_t cds_nt_remainder = 0;     // 0, 1 or 2
    // "ENSP" if the input id matched a protein, "ENST" if it matched a
    // transcript directly. Written to summary so users can audit.
    std::string input_id_type;        // "ENSP" | "ENST" | ""
};

struct UnmappedDomainRow {
    std::string input_id;
    std::string protein_id;
    uint32_t aa_start;
    uint32_t aa_end;
    std::string domain_id;
    std::string reason;
};

struct DomainResult {
    ProteinDomain domain;
    bool mapped;
    bool no_domain_mode;
    UnmappedDomainRow unmapped;

    MappingSummaryRow summary;
    SpanRow span;
    bool has_span;
    std::vector<IsoformSegmentRow> isoform_segments;

    DomainResult() : mapped(false), no_domain_mode(false), has_span(false) {}
};

class DomainMapper {
public:
    DomainMapper(const GTFParser& parser, OutputKind output_kind);

    ErrorCode load_domains(const std::string& bed_filename);
    ErrorCode process_domains();

    // Process queries in chunks of `batch_size`, invoking `on_chunk` once per
    // chunk with the results for that slice of the input. After the callback
    // returns, the chunk is destroyed. Use this path to bound peak RAM at
    // roughly O(batch_size * per-query result size); the one-shot
    // process_domains() keeps O(N) like before.
    //
    // results() is NOT populated by the streaming path; the streaming caller
    // is expected to consume each chunk inside on_chunk. mapped/unmapped
    // counts ARE tracked so the same accessors keep working afterwards.
    ErrorCode process_domains_streaming(
        size_t batch_size,
        const std::function<void(std::vector<DomainResult>&)>& on_chunk);

    const std::vector<DomainResult>& results() const { return results_; }
    size_t domain_count() const { return domains_.size(); }
    size_t mapped_count() const;
    size_t unmapped_count() const { return domain_count() - mapped_count(); }

private:
    const GTFParser& gtf_;
    OutputKind output_kind_;
    std::vector<ProteinDomain> domains_;
    std::vector<DomainResult> results_;
    // Tracked separately so the streaming path (which never populates
    // results_) can still answer mapped_count() / unmapped_count().
    size_t streamed_mapped_ = 0;
    size_t streamed_unmapped_ = 0;
    bool streamed_ = false;

    DomainResult process_one(const ProteinDomain& d) const;
};

#endif // DOMAIN_MAPPER_HPP
