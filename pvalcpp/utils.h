#include <algorithm>
#include <iterator>
#include <limits>
#include <vector>
#include <cassert>

#include "iterator.h"

namespace scoring {

template<class It>
class otherside : public llvm::iterator_facade_base<otherside<It>,
                                                    typename std::iterator_traits<It>::iterator_category,
                                                    typename std::iterator_traits<It>::value_type> {
public:
    using value_type = typename std::iterator_traits<It>::value_type;

    otherside(It it, value_type peptide_mass)
            : it_(it), peptide_mass_{peptide_mass} {}

    value_type operator*() const { return peptide_mass_ - *it_; }
    otherside &operator++() { ++it_; return *this; }
    bool operator==(const otherside &RHS) const { return it_ == RHS.it_; }

    using difference_type = typename std::iterator_traits<It>::difference_type;
    difference_type operator-(const otherside<It> &rhs) const {
      return it_ - rhs.it_;
    }

private:
    value_type peptide_mass_;
    It it_;
};

template <typename It>
otherside<It> make_otherside(It it,
                             typename std::iterator_traits<It>::value_type peptide_mass) {
    return otherside<It>(it, peptide_mass);
}

template<class It>
class charged : public llvm::iterator_facade_base<charged<It>,
                                                  typename std::iterator_traits<It>::iterator_category,
                                                  typename std::iterator_traits<It>::value_type> {
public:
    using value_type = typename std::iterator_traits<It>::value_type;

    charged(It it, unsigned charge = 1)
            : it_(it), inv_charge_{value_type(1) / charge} {}

    value_type operator*() const { return *it_ * inv_charge_; }
    charged &operator++() { ++it_; return *this; }
    bool operator==(const charged &RHS) const { return it_ == RHS.it_; }

    using difference_type = typename std::iterator_traits<It>::difference_type;
    difference_type operator-(const charged<It> &rhs) const {
      return it_ - rhs.it_;
    }

private:
    value_type inv_charge_;
    It it_;
};

template<class It>
charged<It> make_charged(It it, unsigned charge = 1) {
    return charged<It>(it, charge);
}

template<class It>
class sentinelled : public llvm::iterator_facade_base<sentinelled<It>,
                                                      std::forward_iterator_tag,
                                                      typename std::iterator_traits<It>::value_type> {
public:
    using value_type = typename std::iterator_traits<It>::value_type;

    sentinelled(It it, It end,
                value_type sentinel = std::numeric_limits<value_type>::max())
            : it_{it}, end_{end}, infinity_{sentinel} {}

    value_type operator*() const { return it_ != end_ ? *it_ : infinity_; }
    sentinelled &operator++() { if (it_ != end_) ++it_; return *this; }
    sentinelled operator++(int) {
        sentinelled tmp = *this;
        if (it_ != end_) ++it_;
        return tmp;
    }

    bool operator==(const sentinelled &RHS) const { return it_ == RHS.it_; }
private:
    It it_, end_;
    value_type infinity_;
};

template<class It>
sentinelled<It> make_sentinelled(It it, It end) {
    return sentinelled<It>(it, end);
}

template<class It>
std::vector<typename std::iterator_traits<It>::value_type>
add_right_tail_to_spectrum(It spectrum_begin, It spectrum_end,
                           double peptide_mass) {
    using value_type = typename std::iterator_traits<It>::value_type;

    std::vector<value_type> result;
    result.reserve(2 * std::distance(spectrum_begin, spectrum_end));
    result.insert(result.begin(), spectrum_begin, spectrum_end);

    while (spectrum_begin != spectrum_end)
        result.push_back(peptide_mass - *spectrum_begin++);

    std::sort(result.begin(), result.end());
    return result;
}

template<class It, class OutIt>
OutIt add_right_tail_to_spectrum_sorted(It spectrum_begin, It spectrum_end,
                                        OutIt d_out,
                                        double peptide_mass) {
    // assert(std::is_sorted(spectrum_begin, spectrum_end));

    return std::merge(spectrum_begin, spectrum_end,
                      make_otherside(std::reverse_iterator<It>(spectrum_end), peptide_mass),
                      make_otherside(std::reverse_iterator<It>(spectrum_begin), peptide_mass),
                      d_out);
}

template<class It>
std::vector<typename std::iterator_traits<It>::value_type>
add_charged_peaks_to_spectrum(It spectrum_begin, It spectrum_end,
                              unsigned max_charge = 1) {
    using value_type = typename std::iterator_traits<It>::value_type;

    std::vector<value_type> result;
    result.reserve(max_charge * std::distance(spectrum_begin, spectrum_end));
    
    while (spectrum_begin != spectrum_end) {
        double peak = *spectrum_begin++;
        for (size_t charge = 1; charge <= max_charge; ++charge)
            result.push_back(peak / charge);
    }

    std::sort(result.begin(), result.end());

    return result;
}

template<class It, class OutIt>
OutIt add_charged_peaks_to_spectrum_sorted(It spectrum_begin, It spectrum_end,
                                           OutIt d_out,
                                           unsigned max_charge = 1) {
    // assert(std::is_sorted(spectrum_begin, spectrum_end));
    size_t outsz = max_charge * std::distance(spectrum_begin, spectrum_end);
    if (max_charge == 1) {
        return std::copy(spectrum_begin, spectrum_end, d_out);
    } else if (max_charge == 2) {
        auto it1 = make_sentinelled(spectrum_begin, spectrum_end);
        auto it2 = make_sentinelled(make_charged(spectrum_begin, 2), make_charged(spectrum_end, 2));
        for (size_t i = 0; i < outsz; ++i) {
            if (*it1 < *it2) {
                *d_out++ = *it1++;
            } else {
                *d_out++ = *it2++;
            }
        }
    } else if (max_charge == 3) {
        auto it1 = make_sentinelled(spectrum_begin, spectrum_end);
        auto it2 = make_sentinelled(make_charged(spectrum_begin, 2), make_charged(spectrum_end, 2));
        auto it3 = make_sentinelled(make_charged(spectrum_begin, 3), make_charged(spectrum_end, 3));
        for (size_t i = 0; i < outsz; ++i) {
            if (*it1 < *it2 && *it1 < *it3) {
                *d_out++ = *it1++;
            } else if (*it2 < *it3) {
                *d_out++ = *it2++;
            } else {
                *d_out++ = *it3++;
            }
        }
    } else if (max_charge == 4) {
        auto it1 = make_sentinelled(spectrum_begin, spectrum_end);
        auto it2 = make_sentinelled(make_charged(spectrum_begin, 2), make_charged(spectrum_end, 2));
        auto it3 = make_sentinelled(make_charged(spectrum_begin, 3), make_charged(spectrum_end, 3));
        auto it4 = make_sentinelled(make_charged(spectrum_begin, 4), make_charged(spectrum_end, 4));
        for (size_t i = 0; i < outsz; ++i) {
            if (*it1 < *it2 && *it1 < *it3 && *it1 < *it4) {
                *d_out++ = *it1++;
            } else if (*it2 < *it3 && *it2 < *it4) {
                *d_out++ = *it2++;
            } else if (*it3 < *it4) {
                *d_out++ = *it3++;
            } else {
                *d_out++ = *it4++;
            }
        }
    } else {
        std::vector<typename std::iterator_traits<It>::value_type> result;
        while (spectrum_begin != spectrum_end) {
            auto peak = *spectrum_begin++;
            for (size_t charge = 1; charge <= max_charge; ++charge)
                result.push_back(peak / charge);
        }

        std::sort(result.begin(), result.end());
        return std::copy(result.begin(), result.end(), d_out);
    }

    return d_out;
}

} // namespace scoring

    
