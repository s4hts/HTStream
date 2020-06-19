#ifndef MIN_SCREENER_H
#define MIN_SCREENER_H

#include "ioHandler.h"
#include "utils.h"
#include "threadutils.h"
#include "main_template.h"
#include "read.h"
#include "counters.h"
#include "minimizer.h"
#include "phix.h"
#include "../../hts_SeqScreener/src/hts_SeqScreener.h"

#include <algorithm>

extern template class InputReader<SingleEndRead, SingleEndReadFastqImpl>;
extern template class InputReader<PairedEndRead, PairedEndReadFastqImpl>;
extern template class InputReader<PairedEndRead, InterReadImpl>;
extern template class InputReader<ReadBase, TabReadImpl>;

class MinScreener: public MainTemplate<SeqScreenerCounters, MinScreener> {
public:

    MinScreener() {
        program_name = "hts_MinScreener";
        app_description =
            "Adapter Trimmer, trims off adapters by overlapping paired-end reads and\n";
        app_description += "  trimming off overhangs which by definition are adapter sequence in standard\n";
        app_description += "  libraries. SE Reads are trimmed by overlapping the adapter-sequence and trimming off the overlap.";
    }

    void add_extra_options(po::options_description &desc) {
        desc.add_options()
            ("seq,s", po::value<std::string>(), "Please supply a fasta file - default - Phix Sequence - default https://www.ncbi.nlm.nih.gov/nuccore/9626372")
            ("check-read-2,C", po::bool_switch()->default_value(false), "Check R2 as well as R1 (pe)")
            ("kmer,k", po::value<size_t>()->default_value(12)->notifier(boost::bind(&check_range<size_t>, "kmer", _1, 5, 32)), "Kmer size of the lookup table (min 5, max 32)")
            ("window,w", po::value<size_t>()->default_value(12)->notifier(boost::bind(&check_range<size_t>, "kmer", _1, 5, 32)), "Kmer window size, number of consecutive kmers (min 5, max 32)")
            ("percentage-hits,x", po::value<double>()->default_value(.25)->notifier(boost::bind(&check_range<double>, "percentage-hits", _1, 0.0, 1.0)), "Proportion of kmer percentage-hits to sequence need to happen to discard (min 0.0, max 1.0)")
            ("inverse,n", po::bool_switch()->default_value(false), "Output reads that are ABOVE the kmer hit threshold")
            ("record,r", po::bool_switch()->default_value(false), "Only record the reads that pass the kmer hit threshold, output all reads");
    }


    void load_fasta(InputReader<Reference, FastaReadImpl> &faReader) {
        while(faReader.has_next()) {
            auto ref = faReader.next();
            ref_minimizer->find_mins(ReferencePtr(std::move(ref)));
        }

        size_t count = ref_minimizer->get_kmers().bucket_count();
        float load = ref_minimizer->get_kmers().load_factor();
        std::cerr << "reference buckets: " << count << ", elements per bucket: " << load << std::endl;
    }


    size_t check_read(const ReadPtr read) {
        // clear for each read
        read_minimizer->get_kmers().clear();

        read_minimizer->find_mins(read);

        size_t hits = 0;

        for (auto& mini : read_minimizer->get_kmers()) {
            hits += ref_minimizer->get_kmers().count(mini.first);
        }
        // if (hits > 0) {
        //     std::cerr << "found " << hits << " from read " << read->get_seq() << std::endl;
        // }
        return hits;
    }

    template <class T, class Impl>
    void do_app(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, SeqScreenerCounters &counter, const po::variables_map &vm) {

        ref_minimizer = std::unique_ptr<Minimizer<ReferencePtr>>(
            new Minimizer<ReferencePtr>(vm["window"].as<size_t>(), vm["kmer"].as<size_t>()));
        read_minimizer = std::unique_ptr<Minimizer<ReadPtr>>(
            new Minimizer<ReadPtr>(vm["window"].as<size_t>(), vm["kmer"].as<size_t>()));
        WriterHelper writer(pe, se, false);

        std::string lookup_file;
        if (vm.count("seq")) {
            lookup_file = vm["seq"].as<std::string>();
            bi::stream <bi::file_descriptor_source> fa{check_open_r(lookup_file), bi::close_handle};
            InputReader<Reference, FastaReadImpl> faReader(fa);
            load_fasta(faReader);
        } else {
            std::istringstream ins(phixSeqFasta);
            InputReader<Reference, FastaReadImpl> faReader(ins);
            load_fasta(faReader);
        }
        if (ref_minimizer->get_kmers().size() == 0) {
            throw HtsRuntimeException("Fasta reference contains no valid data");
        }


        auto read_visit = make_read_visitor_func(
            [&](SingleEndRead *ser) {
                size_t hits = check_read(ser->get_read_ptr());
                if (hits > 0) {
                    counter.inc_SE_hits();
                }

                counter.output(*ser);
                writer(*ser);
            },
            [&](PairedEndRead *per) {
                size_t hits = check_read(per->get_read_one_ptr());
                hits += check_read(per->get_read_two_ptr());

                if (hits > 0) {
                    counter.inc_PE_hits();
                }

                counter.output(*per);
                writer(*per);
            });

        while(reader.has_next()) {
            auto i = reader.next();
            counter.input(*i);
            i->accept(read_visit);
        }
    }

private:
    std::unique_ptr<Minimizer<ReferencePtr>> ref_minimizer;
    std::unique_ptr<Minimizer<ReadPtr>> read_minimizer;
};
#endif
