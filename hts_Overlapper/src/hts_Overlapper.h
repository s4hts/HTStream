#ifndef OVERLAPPER_H
#define OVERLAPPER_H
//  this is so we can implment hash function for dynamic_bitset
#define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS
#define STARTS 4098

#include "ioHandler.h"
#include <map>
#include <unordered_map>
#include <boost/dynamic_bitset.hpp>
#include <boost/functional/hash.hpp>
#include <algorithm>
#include <bitset>
#include <utility>
#include "utils.h"
#include "main_template.h"
#include "threadutils.h"

typedef std::pair<ReadBasePtr, SingleEndReadPtr> OverlapPair;

class OverlappingCounters : public Counters {
public:
    std::vector<uint_fast64_t> insertLength;

    std::vector<Label> overlapped;

    uint64_t sins = 0;
    uint64_t mins = 0;
    uint64_t lins = 0;

    OverlappingCounters(const std::string &program_name, const po::variables_map &vm) : Counters::Counters(program_name, vm) {

        insertLength.resize(1);

        overlapped.push_back(std::forward_as_tuple("short", sins));
        overlapped.push_back(std::forward_as_tuple("medium", mins));
        overlapped.push_back(std::forward_as_tuple("long", lins));
    }

    using Counters::output;
    void overlap_stats(SingleEndRead &ser, uint_fast64_t origLength) {
        Read &one = ser.non_const_read_one();
        if (one.getLength() < origLength) {
            ++sins; //adapters must be had (short insert)
        } else {
            ++mins; //must be a long insert
        }
        if ( one.getLength() + 1 > insertLength.size() ) {
            insertLength.resize(one.getLength() + 1);
        }
        ++insertLength[one.getLength()];
    }

    void increment_lins()  {
        ++lins;
    }

    void write_out() {

        std::vector<Vector> iLength;
        for (size_t i = 1; i < insertLength.size(); ++i) {
            if (insertLength[i] > 0) {
                iLength.push_back(std::forward_as_tuple(i, insertLength[i]));
            }
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
        start_sublabel("inserts",2);
        write_values(overlapped, 3);
        end_sublabel(2);
        write_vector("overlap_histogram",iLength, 2);
        end_sublabel();

        start_sublabel("Single_end");
        write_values(se, 2);
        end_sublabel();

        start_sublabel("Paired_end");
        write_values(pe, 2);
        start_sublabel("Read1",2);
        write_values(r1, 3);
        end_sublabel(2);
        start_sublabel("Read2",2);
        write_values(r2, 3);
        end_sublabel(2);

        end_sublabel();

        finalize_json();
    }
};

class Overlapper: public MainTemplate<OverlappingCounters, Overlapper> {
public:

    void add_extra_options(po::options_description &desc) {
        setThreadPoolParams(desc);
        // number-of-threads|p

        setDefaultParamsOverlapping(desc);
        // kmer|k ; kmer-offset|r ; max-mismatch-errorDensity|x
        // check-lengths|c ; min-overlap|o

        desc.add_options()
            ("force-pairs,X", po::bool_switch()->default_value(false), "after overlapping a paired end read, split reads in half to output pairs.");
    }

    Overlapper() {
        program_name = "hts_Overlapper";
        app_description =
            "The hts_Overlapper application attempts to overlap paired end reads\n";
        app_description += "  to produce the original transcript, trim adapters, and in some\n";
        app_description += "  cases, correct sequencing errors. single end reads are passed through unchanged.\n";
        app_description += "Reads come in three flavors:\n";
        app_description += "  short: Reads produced from an insert shorter than the longest read\n";
        app_description += "        will result in a single read in the orientation of R1, and have overhanging\n";
        app_description += "        bases (adapters) trimmed to produce a SE read.\n";
        app_description += "  medium: Reads produced from a medium-insert greater than read length, but\n";
        app_description += "        somewhat shorter than 2x read length will produce a SE read in the\n";
        app_description += "        orientation of R1.\n";
        app_description += "  long: Reads produced from long-inserts which do not overlap\n";
        app_description += "  by at least min overlap , resulting in a PE read.\n";
    }

/* Within the overlap if they are the same bp, then add q scores
 * If they are different bp, subtract q scores and take the larger quality bp*/
    SingleEndReadPtr checkIfOverlap(Read &r1, Read &r2, size_t loc1, size_t loc2, const double misDensity, const size_t &mismatch, const size_t minOverlap) {
        size_t minLoc = std::min(loc1, loc2);
        size_t loc1_t = loc1 - minLoc;
        size_t loc2_t = loc2 - minLoc;
        size_t r1_len = r1.getLength();
        size_t r2_len = r2.getLength();

        size_t maxLoop = std::min(r1_len - loc1_t, r2_len - loc2_t);
        size_t maxMis = std::min(mismatch, static_cast<size_t>(maxLoop * misDensity));

        const std::string &seq1 = r1.get_seq();
        const std::string &seq2 = r2.get_seq_rc();

        auto i1 = seq1.begin();
        std::advance(i1, loc1_t);
        auto i2 = seq2.begin();
        std::advance(i2, loc2_t);
        if (maxLoop < minOverlap || !threshold_mismatches( i1, i2 , maxLoop, maxMis ) ) {
            return nullptr;
        }

        size_t read1_bp;
        size_t read2_bp;

        const std::string &qual1 = r1.get_qual();
        const std::string &qual2 = r2.get_qual_rc();

        std::string finalSeq;
        std::string finalQual;

        char bp;
        char qual;

        for (size_t i = 0; i < maxLoop; ++i) {
            read1_bp = loc1_t + i;
            read2_bp = loc2_t + i;

            if (seq1[read1_bp] == seq2[read2_bp]) {
                bp = seq1[read1_bp];
                qual = static_cast<char>(std::min(qual1[read1_bp] + qual2[read2_bp] - 33, 40 + 33)); //addition of qual (minus just one of the ascii values
            } else {
                bp = qual1[read1_bp] >= qual2[read2_bp] ? seq1[read1_bp] : seq2[read2_bp];
                qual = static_cast<char>(std::max(qual1[read1_bp] - qual2[read2_bp] + 33, 1 + 33));
            }
            finalSeq += bp;
            finalQual += qual;
        }

        /*R2 is always shorter or equal to R1*/

        if ( loc1_t > loc2_t ) { //We are going to snag the start if it is a lin (notice, if they are both zero then we dont' want to snag anything
            finalSeq = seq1.substr(0, loc1_t) + finalSeq;
            finalQual = qual1.substr(0, loc1_t) + finalQual;
        }

        if (r1_len - loc1_t < r2_len - loc2_t) {
            finalSeq += seq2.substr(maxLoop, r2_len - maxLoop);
            finalQual += qual2.substr(maxLoop, r2_len - maxLoop);
        }

        Read overlap(finalSeq, finalQual, r1.get_id_tab("1"));
        overlap.join_comment(r1.get_comment());
        overlap.join_comment(r2.get_comment());
        SingleEndReadPtr overlap_se(new SingleEndRead(overlap));
        return overlap_se;
    }

/*Because of the way overlapping works, you only need to check the ends of the shorter read*/
    SingleEndReadPtr getOverlappedReads(Read &r1, Read &r2, const seqLookup &seq1Map,  const double misDensity, const size_t mismatch, const size_t &minOver, const size_t &checkLengths, const size_t kmer) {
        std::string seq2 = r2.get_seq_rc();
        for (size_t bp = 0; bp < (checkLengths - kmer); ++bp) {
            /*Do a quick check if the shorter read kmer shows up in longer read (read 2)
             * If it does, then try the brute force approach*/
            auto test = seq1Map.equal_range(seq2.substr(bp, kmer));
            for (auto it = test.first; it != test.second; ++it) {
                SingleEndReadPtr overlapped = checkIfOverlap(r1, r2, it->second, bp, misDensity, mismatch, minOver);
                if (overlapped != nullptr) {
                    return overlapped;
                }
            }
        }
        for (size_t bp = seq2.length() - (checkLengths + kmer); bp <= (seq2.length() - kmer) ; ++bp) {
            /*Do a quick check if the shorter read kmer shows up in longer read (read 2)
             * If it does, then try the brute force approach*/
            auto test = seq1Map.equal_range(seq2.substr(bp, kmer));
            for (auto it = test.first; it != test.second; ++it) {
                SingleEndReadPtr overlapped = checkIfOverlap(r1, r2, it->second, bp, misDensity, mismatch, minOver);
                if (overlapped != nullptr) {
                    return overlapped;
                }
            }
        }
        return nullptr;
    }

    OverlapPair check_read(PairedEndReadPtr pe, const double misDensity, const size_t &mismatch, const size_t &minOver, const size_t &checkLengths, const size_t kmer, const size_t kmerOffset) {

        Read &r1 = pe->non_const_read_one();
        Read &r2 = pe->non_const_read_two();

        bool swapped = false;
        /* Read1 needs to always be longer than Read 2 */
        if (r1.getLength() < r2.getLength()) {
            std::swap(r1, r2);
            swapped = true;
        }
        /* checkL needs to be as long as or longer than the shortest read */
        size_t checkL = std::min(r2.getLength(), checkLengths);
        /* checkL needs to be as long as or longer than the shortest read */
        size_t kkmer = std::min(r2.getLength(), kmer);
        /* Create a map with non-overlapping kmers */
        seqLookup mOne = readOneMap(r1.get_seq(), kkmer, kmerOffset);
        /* returns null if no much
         * r1 and r2 and passed by ref in case only adapter trimming is on */
        SingleEndReadPtr overlapped = getOverlappedReads(r1, r2, mOne, misDensity, mismatch, minOver, checkL, kkmer) ;
        if (!overlapped && swapped){
            std::swap(r1, r2);
        } else if (swapped){
            overlapped->set_read_rc();
        }

        return std::make_pair(std::static_pointer_cast<ReadBase>(pe), overlapped);
    }

    void writer_thread(std::shared_ptr<OutputWriter> pe,  std::shared_ptr<OutputWriter> se, OverlappingCounters &counter, threadsafe_queue<std::future<OverlapPair>> &futures, bool forcePair) {

        WriterHelper writer(pe, se);

        SingleEndReadPtr overlapped;

        auto read_visit = make_read_visitor_func(
            [&](SingleEndRead *ser) {
                counter.output(*ser);
                writer(*ser);
            },
            [&](PairedEndRead *per) {

                if (!overlapped) {
                    counter.increment_lins();
                    counter.output(*per);
                    writer(*per); //write out as is
                } else if (overlapped) { //if there is an overlap
                    unsigned int origLength = std::max(unsigned(per->non_const_read_one().getLength()),unsigned(per->non_const_read_two().getLength()));
                    counter.overlap_stats(*overlapped, origLength);
                    if (forcePair){
                        Read overlappedRead = overlapped->non_const_read_one();
                        Read &or1 = per->non_const_read_one();
                        Read &or2 = per->non_const_read_two();
                        double mid = ceil(overlappedRead.getLengthTrue() / 2.0);
                        Read r1(overlappedRead.get_sub_seq().substr(0,mid),overlappedRead.get_sub_qual().substr(0,mid),overlappedRead.get_id_tab("1"));
                        r1.join_comment(or1.get_comment());
                        Read r2(overlappedRead.get_sub_seq().substr(mid),overlappedRead.get_sub_qual().substr(mid),overlappedRead.get_id_tab("1"));
                        r2.join_comment(or2.get_comment());
                        PairedEndRead newper(r1, r2);
                        counter.output(newper);
                        writer(newper); //write out as is
                    } else {
                        counter.output(*overlapped);
                        writer(*overlapped);
                    }
                }
            });


        while(!futures.is_done()) {
            std::future<OverlapPair> fread;
            futures.wait_and_pop(fread);

            OverlapPair opair = fread.get();

            // null read indicates all done
            if (!opair.first) {
                futures.set_done();
                return;
            }
            overlapped = opair.second;
            opair.first->accept(read_visit);
        }
    }

/*This is the helper class for overlap
 * The idea is in the wet lab, they set up sequences to sequences toward each other
 * Sometimes, these squences can overlap
 * There are two cases in which they can overlap
 * They can either overlap just barely on the ends we call - these mins (medium insert)
 * They can also overlap way to much to the point they have adapters
 * in the read - or a sin (short insert) [adapter are overhangs in an overlap]
 * With a min it is useful to have a higher confidence in the bases in the overlap and longer read
 * With a sin it is useful to have the higher confidence as well as removing the adapters*/
    template <class T, class Impl>
    void do_app(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, OverlappingCounters &counter, const po::variables_map &vm) {
        const double misDensity = vm["max-mismatch-errorDensity"].as<double>();
        const size_t mismatch = vm["max-mismatch"].as<size_t>();
        const size_t minOver = vm["min-overlap"].as<size_t>();
        const size_t checkLengths = vm["check-lengths"].as<size_t>();
        const size_t kmer = vm["kmer"].as<size_t>();
        const size_t kmerOffset = vm["kmer-offset"].as<size_t>();
        bool forcePair = vm["force-pairs"].as<bool>();

        size_t num_threads = vm["number-of-threads"].as<size_t>();

        threadsafe_queue<std::future<OverlapPair>> futures(50000);
        thread_pool threads(50000, num_threads);

        std::thread output_thread([=, &counter, &futures]() mutable { writer_thread(pe, se, counter, futures, forcePair); });

        auto read_visit = make_read_visitor_func(
            [&](SingleEndRead *ser) {
                SingleEndReadPtr sser(ser);
                futures.push(threads.submit([=]() { return std::make_pair(std::static_pointer_cast<ReadBase>(sser), SingleEndReadPtr()); }));
            },
            [&](PairedEndRead *per) {
                futures.push(threads.submit([=]() mutable {
                                                return check_read(PairedEndReadPtr(per), misDensity, mismatch, minOver, checkLengths, kmer, kmerOffset); }));
            });
        // thread_guard must be declared last
        thread_guard tg(output_thread);

        while(reader.has_next()) {
            auto i = reader.next();
            counter.input(*i);

            // release so we convert to shared_ptr later
            auto p = i.release();
            p->accept(read_visit);
        }
        // null ptr indicates end of processing
        futures.push(threads.submit([]() { return std::make_pair(ReadBasePtr(), SingleEndReadPtr()); }));
    }
};
#endif
