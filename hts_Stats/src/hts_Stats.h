#ifndef STATS_H
#define STATS_H
//  this is so we can implment hash function for dynamic_bitset
#define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS

#include "ioHandler.h"
#include "utils.h"
#include "main_template.h"
#include "threadutils.h"

#include <array>
#include <deque>
#include <map>
#include <unordered_map>
#include <algorithm>

extern template class InputReader<SingleEndRead, SingleEndReadFastqImpl>;
extern template class InputReader<PairedEndRead, PairedEndReadFastqImpl>;
extern template class InputReader<PairedEndRead, InterReadImpl>;
extern template class InputReader<ReadBase, TabReadImpl>;

class StatsCounters : public Counters {

public:
    typedef std::array<uint_fast64_t, 5> BaseCycle;
    typedef std::array<uint_fast64_t, QUAL_MAX> QualityCycle;

    Vec R1_Length;
    Vec R2_Length;
    Vec SE_Length;

    std::vector<BaseCycle> R1_bases;
    std::vector<BaseCycle> R2_bases;
    std::vector<BaseCycle> SE_bases;

    std::vector<QualityCycle> R1_qualities;
    std::vector<QualityCycle> R2_qualities;
    std::vector<QualityCycle> SE_qualities;

    std::vector<Label> bases;

    std::array<uint64_t, 5> base_counts{{0, 0, 0, 0, 0}};
    uint64_t &A;
    uint64_t &C;
    uint64_t &G;
    uint64_t &T;
    uint64_t &N;

    uint64_t SE_bQ30 = 0;

    uint64_t R1_bQ30 = 0;
    uint64_t R2_bQ30 = 0;

    size_t qual_offset;
    size_t max_quality_seen = 0;
    std::array<unsigned char, 256> base_lookup;

    StatsCounters(const std::string &program_name, const po::variables_map &vm) :
        Counters::Counters(program_name, vm),
        A(base_counts[0]),
        C(base_counts[1]),
        G(base_counts[2]),
        T(base_counts[3]),
        N(base_counts[4]) {
        R1_Length.resize(1);
        R2_Length.resize(1);
        SE_Length.resize(1);
        qual_offset = vm.count("qual-offset") ? vm["qual-offset"].as<size_t>() : DEFAULT_QUAL_OFFSET;
        base_lookup.fill(255);
        base_lookup[static_cast<unsigned char>('A')] = 0;
        base_lookup[static_cast<unsigned char>('C')] = 1;
        base_lookup[static_cast<unsigned char>('G')] = 2;
        base_lookup[static_cast<unsigned char>('T')] = 3;
        base_lookup[static_cast<unsigned char>('N')] = 4;
        se.push_back(std::forward_as_tuple("total_Q30_basepairs", SE_bQ30));
        r1.push_back(std::forward_as_tuple("total_Q30_basepairs", R1_bQ30));
        r2.push_back(std::forward_as_tuple("total_Q30_basepairs", R2_bQ30));

        bases.push_back(std::forward_as_tuple("A", A));
        bases.push_back(std::forward_as_tuple("C", C));
        bases.push_back(std::forward_as_tuple("G", G));
        bases.push_back(std::forward_as_tuple("T", T));
        bases.push_back(std::forward_as_tuple("N", N));
    }
    virtual ~StatsCounters() {}

    static void merge_lengths(Vec &target, const Vec &source) {
        if (target.size() < source.size()) {
            target.resize(source.size());
        }
        for (size_t i = 0; i < source.size(); ++i) {
            target[i] += source[i];
        }
    }

    template <size_t N>
    static void merge_cycles(std::vector<std::array<uint_fast64_t, N> > &target, const std::vector<std::array<uint_fast64_t, N> > &source) {
        if (target.size() < source.size()) {
            const size_t old_size = target.size();
            target.resize(source.size());
            for (size_t i = old_size; i < target.size(); ++i) {
                target[i].fill(0);
            }
        }
        for (size_t cycle = 0; cycle < source.size(); ++cycle) {
            for (size_t row = 0; row < N; ++row) {
                target[cycle][row] += source[cycle][row];
            }
        }
    }

    void merge_from(const StatsCounters &other) {
        TotalFragmentsInput += other.TotalFragmentsInput;
        TotalFragmentsOutput += other.TotalFragmentsOutput;
        TotalBasepairsInput += other.TotalBasepairsInput;
        TotalBasepairsOutput += other.TotalBasepairsOutput;

        SE_In += other.SE_In;
        SE_Out += other.SE_Out;
        SE_BpLen_In += other.SE_BpLen_In;
        SE_BpLen_Out += other.SE_BpLen_Out;

        PE_In += other.PE_In;
        PE_Out += other.PE_Out;
        R1_BpLen_In += other.R1_BpLen_In;
        R1_BpLen_Out += other.R1_BpLen_Out;
        R2_BpLen_In += other.R2_BpLen_In;
        R2_BpLen_Out += other.R2_BpLen_Out;

        merge_lengths(R1_Length, other.R1_Length);
        merge_lengths(R2_Length, other.R2_Length);
        merge_lengths(SE_Length, other.SE_Length);

        merge_cycles(R1_bases, other.R1_bases);
        merge_cycles(R2_bases, other.R2_bases);
        merge_cycles(SE_bases, other.SE_bases);
        merge_cycles(R1_qualities, other.R1_qualities);
        merge_cycles(R2_qualities, other.R2_qualities);
        merge_cycles(SE_qualities, other.SE_qualities);

        for (size_t i = 0; i < base_counts.size(); ++i) {
            base_counts[i] += other.base_counts[i];
        }
        SE_bQ30 += other.SE_bQ30;
        R1_bQ30 += other.R1_bQ30;
        R2_bQ30 += other.R2_bQ30;
        max_quality_seen = std::max(max_quality_seen, other.max_quality_seen);
    }

    template <size_t N>
    static Mat cycles_to_mat(const std::vector<std::array<uint_fast64_t, N> >& cycles) {
        Mat out;
        out.reserve(cycles.size());
        for (const auto& cycle : cycles) {
            out.push_back(Vec(cycle.begin(), cycle.end()));
        }
        return out;
    }

    void read_stats(Read &r, Vec &Length, std::vector<BaseCycle> &read_bases, std::vector<QualityCycle> &read_qualities, uint64_t &read_bQ30) {
        read_stats_strings(r.get_seq(), r.get_qual(), Length, read_bases, read_qualities, read_bQ30);
    }

    void read_stats_strings(const std::string &seq, const std::string &qual, Vec &Length, std::vector<BaseCycle> &read_bases, std::vector<QualityCycle> &read_qualities, uint64_t &read_bQ30) {
        const size_t length = seq.size();
        // Size histogram per read
        if ( length + 1 > Length.size() ) {
            Length.resize(length + 1);
        }
        ++Length[length];
        // READ Base and Quality stats
        // update size of base and Q score matrix if needed
        while(read_bases.size() < length) {
            read_bases.emplace_back();
            read_bases.back().fill(0);
            read_qualities.emplace_back();
            read_qualities.back().fill(0);
        }
        uint64_t q30bases=0;
        for (size_t index = 0; index < length; ++index) {
            // bases
            const unsigned char base_index = base_lookup[static_cast<unsigned char>(seq[index])];
            if (base_index == 255) {
                throw HtsRuntimeException(std::string("Unknown bp in stats counter: ") + seq[index]);
            }
            ++base_counts[base_index];
            ++read_bases[index][base_index];
            // qualities
            const unsigned char qscore = static_cast<unsigned char>(qual[index]);
            if (qscore >= qual_offset) {
                const size_t qscore_int = qscore - qual_offset;
                if (qscore_int < QUAL_MAX) {
                    if (qscore_int > max_quality_seen) {
                        max_quality_seen = qscore_int;
                    }
                    q30bases += qscore_int >= 30;
                    ++read_qualities[index][qscore_int];
                }
            }
        }
        read_bQ30 += q30bases;
    }

    void count_paired_fastq_record(const std::string &seq1, const std::string &qual1, const std::string &seq2, const std::string &qual2) {
        ++PE_In;
        ++PE_Out;
        ++TotalFragmentsInput;
        ++TotalFragmentsOutput;

        R1_BpLen_In += seq1.size();
        R1_BpLen_Out += seq1.size();
        R2_BpLen_In += seq2.size();
        R2_BpLen_Out += seq2.size();
        TotalBasepairsInput += seq1.size() + seq2.size();
        TotalBasepairsOutput += seq1.size() + seq2.size();

        read_stats_strings(seq1, qual1, R1_Length, R1_bases, R1_qualities, R1_bQ30);
        read_stats_strings(seq2, qual2, R2_Length, R2_bases, R2_qualities, R2_bQ30);
    }

    using Counters::output;
    void output(PairedEndRead &per) {
        Counters::output(per);
        read_stats(per.non_const_read_one(), R1_Length, R1_bases, R1_qualities, R1_bQ30);
        read_stats(per.non_const_read_two(), R2_Length, R2_bases, R2_qualities, R2_bQ30);
    }

    void output(SingleEndRead &ser) {
        Counters::output(ser);
        read_stats(ser.non_const_read_one(), SE_Length, SE_bases, SE_qualities, SE_bQ30);
    }

    void write_out() {
        std::vector<Vector> iSE_Length;
        for (size_t i = 0; i < SE_Length.size(); ++i) {
            if (SE_Length[i] > 0) {
                iSE_Length.push_back(std::forward_as_tuple(i, SE_Length[i]));
            }
        }

        std::vector<Vector> iR1_Length;
        for (size_t i = 0; i < R1_Length.size(); ++i) {
            if (R1_Length[i] > 0) {
                iR1_Length.push_back(std::forward_as_tuple(i, R1_Length[i]));
            }
        }

        std::vector<Vector> iR2_Length;
        for (size_t i = 0; i < R2_Length.size(); ++i) {
            if (R2_Length[i] > 0) {
                iR2_Length.push_back(std::forward_as_tuple(i, R2_Length[i]));
            }
        }

        std::vector<std::string> ind_se;
        for (size_t j = 1; j <= SE_bases.size(); j++){
          ind_se.push_back(std::to_string((int)j));
        }
        std::vector<std::string> ind_pe1;
        for (size_t j = 1; j <= R1_bases.size(); j++){
          ind_pe1.push_back(std::to_string((int)j));
        }
        std::vector<std::string> ind_pe2;
        for (size_t j = 1; j <= R2_bases.size(); j++){
          ind_pe2.push_back(std::to_string((int)j));
        }
        std::vector<std::string> b{ "A", "C", "G", "T", "N"};
        std::vector<std::string> q;
        for (size_t j = 0; j <= max_quality_seen && j < QUAL_MAX; j++){
          q.push_back(std::to_string((int)j));
        }

        initialize_json();

        start_sublabel("Program_details");
        write_values(pd, 2);
        start_sublabel("options",2);
        write_options(3);
        end_sublabel(2);
        end_sublabel();

        start_sublabel("Fragment");
        write_values(fragment, 2);
        start_sublabel("base_composition", 2);
        write_values(bases, 3);
        end_sublabel(2);
        end_sublabel();

        start_sublabel("Single_end");
        write_values(se, 2);
        write_vector("readlength_histogram",iSE_Length, 2);
        write_matrix("base_by_cycle", cycles_to_mat(SE_bases), b, ind_se, 0, 2);
        write_matrix("qualities_by_cycle", cycles_to_mat(SE_qualities), q, ind_se, 0, 2);
        end_sublabel();

        start_sublabel("Paired_end");
        write_values(pe, 2);
        start_sublabel("Read1",2);
        write_values(r1, 3);
        write_vector("readlength_histogram",iR1_Length, 3);
        write_matrix("base_by_cycle", cycles_to_mat(R1_bases), b, ind_pe1, 0, 3);
        write_matrix("qualities_by_cycle", cycles_to_mat(R1_qualities), q, ind_pe1, 0, 3);
        end_sublabel(2);
        start_sublabel("Read2",2);
        write_values(r2, 3);
        write_vector("readlength_histogram",iR2_Length, 3);
        write_matrix("base_by_cycle", cycles_to_mat(R2_bases), b, ind_pe2, 0, 3);
        write_matrix("qualities_by_cycle", cycles_to_mat(R2_qualities), q, ind_pe2, 0, 3);
        end_sublabel(2);
        end_sublabel();

        finalize_json();
    }

};

class Stats: public MainTemplate<StatsCounters, Stats> {
public:

    struct FastqPairRecord {
        std::string seq1;
        std::string qual1;
        std::string seq2;
        std::string qual2;
    };

    Stats() {
        program_name = "hts_Stats";
        app_description =
            "The hts_Stats app produce basic statistics about the reads in a dataset.\n";
        app_description += "  Including the basepair composition and number of bases Q30.";
    }

    void add_extra_options(po::options_description &desc) {
        setThreadPoolParams(desc);
        desc.add_options()
            ("no-output", po::bool_switch()->default_value(false), "Only write stats JSON; do not pass reads through to stdout or output files");
    }

    template <class T>
    struct BatchStatsResult {
        std::shared_ptr<StatsCounters> counters;
        std::shared_ptr<std::vector<std::unique_ptr<T> > > batch;
    };

    template <class T>
    std::shared_ptr<BatchStatsResult<T> > process_batch(std::shared_ptr<std::vector<std::unique_ptr<T> > > batch, size_t qual_offset) {
        po::variables_map empty_vm;
        std::shared_ptr<StatsCounters> local(new StatsCounters(program_name, empty_vm));
        local->qual_offset = qual_offset;
        for (auto &read : *batch) {
            local->input(*read);
            local->output(*read);
        }
        std::shared_ptr<BatchStatsResult<T> > result(new BatchStatsResult<T>());
        result->counters = local;
        result->batch = batch;
        return result;
    }

    static bool load_fastq_record(std::istream &input, std::string &seq, std::string &qual) {
        std::string id;
        while(std::getline(input, id) && id.size() < 1) {
        }
        if (!input && id.size() < 1) {
            return false;
        }
        if (id.size() < 1 || id[0] != '@') {
            throw HtsIOException("id line did not begin with @");
        }
        std::string id2;
        if (!std::getline(input, seq) || !std::getline(input, id2) || !std::getline(input, qual)) {
            throw HtsIOException("incomplete FASTQ record");
        }
        if (id2.size() < 1 || id2[0] != '+') {
            throw HtsIOException("invalid id2 line did not begin with +");
        }
        if (qual.size() != seq.size()) {
            throw HtsIOException("qual string not the same length as sequence");
        }
        return true;
    }

    std::shared_ptr<StatsCounters> process_fastq_pair_batch(std::shared_ptr<std::vector<FastqPairRecord> > batch, size_t qual_offset) {
        po::variables_map empty_vm;
        std::shared_ptr<StatsCounters> local(new StatsCounters(program_name, empty_vm));
        local->qual_offset = qual_offset;
        for (const auto &record : *batch) {
            local->count_paired_fastq_record(record.seq1, record.qual1, record.seq2, record.qual2);
        }
        return local;
    }

    void submit_fastq_pair_batch(std::deque<std::future<std::shared_ptr<StatsCounters> > > &futures,
                                 thread_pool &threads,
                                 StatsCounters& counters,
                                 std::shared_ptr<std::vector<FastqPairRecord> > batch,
                                 size_t max_pending,
                                 size_t qual_offset) {
        futures.push_back(threads.submit([this, batch, qual_offset]() {
            return process_fastq_pair_batch(batch, qual_offset);
        }));
        while (futures.size() >= max_pending) {
            counters.merge_from(*futures.front().get());
            futures.pop_front();
        }
    }

    bool do_read1_read2_files(const std::vector<std::string> &read1_files, const std::vector<std::string> &read2_files,
                              std::shared_ptr<OutputWriter>, std::shared_ptr<OutputWriter>,
                              StatsCounters& counters, const po::variables_map &vm) {
        if (!vm["no-output"].as<bool>()) {
            return false;
        }

        const size_t num_threads = vm["number-of-threads"].as<size_t>();
        const size_t batch_size = 8192;
        const size_t qual_offset = counters.qual_offset;

        if (num_threads <= 1) {
            for (size_t i = 0; i < read1_files.size(); ++i) {
                bi::stream<bi::file_descriptor_source> is1{check_open_r(read1_files[i]), bi::close_handle};
                bi::stream<bi::file_descriptor_source> is2{check_open_r(read2_files[i]), bi::close_handle};
                std::string seq1, qual1, seq2, qual2;
                while (load_fastq_record(is1, seq1, qual1)) {
                    if (!load_fastq_record(is2, seq2, qual2)) {
                        throw HtsIOException("read2 input ended before read1 input");
                    }
                    counters.count_paired_fastq_record(seq1, qual1, seq2, qual2);
                }
                if (load_fastq_record(is2, seq2, qual2)) {
                    throw HtsIOException("read1 input ended before read2 input");
                }
            }
            return true;
        }

        const size_t max_pending = std::max<size_t>(num_threads * 4, 1);
        thread_pool threads(max_pending, num_threads);
        std::deque<std::future<std::shared_ptr<StatsCounters> > > futures;

        for (size_t i = 0; i < read1_files.size(); ++i) {
            bi::stream<bi::file_descriptor_source> is1{check_open_r(read1_files[i]), bi::close_handle};
            bi::stream<bi::file_descriptor_source> is2{check_open_r(read2_files[i]), bi::close_handle};
            std::shared_ptr<std::vector<FastqPairRecord> > batch(new std::vector<FastqPairRecord>());
            batch->reserve(batch_size);
            std::string seq1, qual1, seq2, qual2;
            while (load_fastq_record(is1, seq1, qual1)) {
                if (!load_fastq_record(is2, seq2, qual2)) {
                    throw HtsIOException("read2 input ended before read1 input");
                }
                FastqPairRecord record;
                record.seq1.swap(seq1);
                record.qual1.swap(qual1);
                record.seq2.swap(seq2);
                record.qual2.swap(qual2);
                batch->push_back(std::move(record));
                if (batch->size() == batch_size) {
                    submit_fastq_pair_batch(futures, threads, counters, batch, max_pending, qual_offset);
                    batch.reset(new std::vector<FastqPairRecord>());
                    batch->reserve(batch_size);
                }
            }
            if (load_fastq_record(is2, seq2, qual2)) {
                throw HtsIOException("read1 input ended before read2 input");
            }
            if (!batch->empty()) {
                submit_fastq_pair_batch(futures, threads, counters, batch, max_pending, qual_offset);
            }
        }

        while (!futures.empty()) {
            counters.merge_from(*futures.front().get());
            futures.pop_front();
        }
        return true;
    }

    template <class T, class Impl>
    void do_parallel_stats(InputReader<T, Impl> &reader, WriterHelper *writer, StatsCounters& counters, const po::variables_map &vm) {
        const size_t num_threads = vm["number-of-threads"].as<size_t>();
        const size_t batch_size = 4096;
        const size_t max_pending = std::max<size_t>(num_threads * 4, 1);
        const size_t qual_offset = counters.qual_offset;
        thread_pool threads(max_pending, num_threads);
        std::deque<std::future<std::shared_ptr<BatchStatsResult<T> > > > futures;

        auto finish_front = [&futures, &counters, writer]() {
            std::shared_ptr<BatchStatsResult<T> > result = futures.front().get();
            futures.pop_front();
            counters.merge_from(*result->counters);
            if (writer) {
                for (auto &read : *result->batch) {
                    (*writer)(*read);
                }
            }
        };

        auto submit_batch = [this, &futures, &threads, &finish_front, max_pending, qual_offset](std::shared_ptr<std::vector<std::unique_ptr<T> > > batch) {
            futures.push_back(threads.submit([this, batch, qual_offset]() {
                return process_batch<T>(batch, qual_offset);
            }));
            while (futures.size() >= max_pending) {
                finish_front();
            }
        };

        std::shared_ptr<std::vector<std::unique_ptr<T> > > batch(new std::vector<std::unique_ptr<T> >());
        batch->reserve(batch_size);
        while(reader.has_next()) {
            batch->push_back(reader.next());
            if (batch->size() == batch_size) {
                submit_batch(batch);
                batch.reset(new std::vector<std::unique_ptr<T> >());
                batch->reserve(batch_size);
            }
        }
        if (!batch->empty()) {
            submit_batch(batch);
        }
        while (!futures.empty()) {
            finish_front();
        }
    }

    template <class T, class Impl>
    void do_app(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, StatsCounters& counters, const po::variables_map &vm) {

        const bool no_output = vm["no-output"].as<bool>();
        WriterHelper writer(pe, se, false);
        if (vm["number-of-threads"].as<size_t>() > 1) {
            do_parallel_stats(reader, no_output ? nullptr : &writer, counters, vm);
            return;
        }

        while(reader.has_next()) {
            auto i = reader.next();
            counters.input(*i);
            counters.output(*i);
            if (!no_output) {
                writer(*i);
            }
        }
    }
};
#endif
