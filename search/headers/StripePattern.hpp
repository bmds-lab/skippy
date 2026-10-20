#ifndef STRIPE_PATTERN_HPP
#define STRIPE_PATTERN_HPP

#include <string>
#include <vector>

// Forward packing - works for any weight (16 or 24+)
inline bool packStripeKmer(std::string_view window, uint64_t& packed)
{
    if ((int)window.size() < gPattern.window)
        return false;

    if (gPattern.use64)
        return gPattern.encode64(window.data(), packed);
    else {
        uint32_t p32 = 0;
        bool ok = gPattern.encode32(window.data(), p32);
        packed = p32;
        return ok;
    }
}

// RC packing
inline bool packStripeKmerRC(const char* wptr, uint64_t& packed)
{
    if (gPattern.use64)
        return gPattern.encode_rc64(wptr, packed);
    else {
        uint32_t p32 = 0;
        bool ok = gPattern.encode_rc32(wptr, p32);
        packed = p32;
        return ok;
    }
}

class StripePattern {
public:
    StripePattern();

    bool parse(const std::string& pattern);

    const std::string& patternString() const;

    const std::vector<int>& blockSizes() const;

    const std::vector<int>& informativeOffsets() const;

    int informativeBases() const;

    int span() const;

    void print() const;

private:
    std::string patternString_;

    std::vector<int> blockSizes_;

    std::vector<int> informativeOffsets_;

    int informativeBases_;

    int span_;
};

#endif
