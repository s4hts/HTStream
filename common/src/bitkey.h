#ifndef BITKEY_H
#define BITKEY_H

#include <algorithm>
#include <cstdint>
#include <ostream>
#include <string>
#include <vector>
#include <boost/optional.hpp>

class BitKey {
public:
    BitKey() : base_count(0) {}

    size_t size() const { return base_count * 2; }
    size_t bases() const { return base_count; }

    bool operator==(const BitKey& other) const {
        return base_count == other.base_count && blocks == other.blocks;
    }

    bool operator!=(const BitKey& other) const {
        return !(*this == other);
    }

    bool operator<(const BitKey& other) const {
        if (base_count != other.base_count) {
            return base_count < other.base_count;
        }
        return blocks < other.blocks;
    }

    bool operator>(const BitKey& other) const {
        return other < *this;
    }

    void append_code(uint8_t code) {
        const size_t bit_offset = (base_count * 2) % BITS_PER_BLOCK;
        if (bit_offset == 0) {
            blocks.push_back(0);
        }
        blocks.back() |= (static_cast<uint64_t>(code) << (BITS_PER_BLOCK - bit_offset - 2));
        ++base_count;
    }

    void append(const BitKey& other) {
        for (size_t i = 0; i < other.base_count; ++i) {
            append_code(other.code_at(i));
        }
    }

    uint8_t code_at(size_t idx) const {
        const size_t bit_pos = idx * 2;
        const size_t block_idx = bit_pos / BITS_PER_BLOCK;
        const size_t bit_offset = bit_pos % BITS_PER_BLOCK;
        return static_cast<uint8_t>((blocks[block_idx] >> (BITS_PER_BLOCK - bit_offset - 2)) & 0x3);
    }

    std::string to_string() const {
        std::string out;
        out.resize(base_count);
        for (size_t i = 0; i < base_count; ++i) {
            switch (code_at(i)) {
            case 0:
                out[i] = 'A';
                break;
            case 1:
                out[i] = 'C';
                break;
            case 2:
                out[i] = 'G';
                break;
            default:
                out[i] = 'T';
                break;
            }
        }
        return out;
    }

    BitKey reverse_complement() const {
        BitKey out;
        for (size_t i = 0; i < base_count; ++i) {
            out.append_code(code_at(base_count - i - 1) ^ 0x3);
        }
        return out;
    }

    static boost::optional<BitKey> from_sequence(const std::string& seq) {
        BitKey out;
        for (const char c : seq) {
            switch (c) {
            case 'A':
            case 'a':
                out.append_code(0);
                break;
            case 'C':
            case 'c':
                out.append_code(1);
                break;
            case 'G':
            case 'g':
                out.append_code(2);
                break;
            case 'T':
            case 't':
                out.append_code(3);
                break;
            default:
                return boost::none;
            }
        }
        return out;
    }

    static boost::optional<BitKey> from_subsequence(const std::string& seq, size_t start, size_t length) {
        BitKey out;
        const size_t end = start + length;
        for (size_t i = start; i < end; ++i) {
            switch (seq[i]) {
            case 'A':
            case 'a':
                out.append_code(0);
                break;
            case 'C':
            case 'c':
                out.append_code(1);
                break;
            case 'G':
            case 'g':
                out.append_code(2);
                break;
            case 'T':
            case 't':
                out.append_code(3);
                break;
            default:
                return boost::none;
            }
        }
        return out;
    }

    friend struct BitKeyHash;

private:
    static const size_t BITS_PER_BLOCK = 64;
    size_t base_count;
    std::vector<uint64_t> blocks;
};

struct BitKeyHash {
    std::size_t operator()(const BitKey& key) const {
        std::size_t seed = key.base_count;
        for (uint64_t block : key.blocks) {
            seed ^= std::hash<uint64_t>()(block) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

inline std::ostream& operator<<(std::ostream& out, const BitKey& key) {
    out << key.to_string();
    return out;
}

#endif
