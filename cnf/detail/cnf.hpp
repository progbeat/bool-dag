#ifndef CNF_DETAIL_CNF_HPP_
#define CNF_DETAIL_CNF_HPP_

#include <unordered_set>

namespace cnf {
namespace detail {

bool has_conflict(const std::unordered_set<int> &vars)
{
	for (auto i = vars.begin(); i != vars.end(); ++i)
		if (*i < 0 && vars.count(~*i))
			return true;
	return false;
}

} // namespace detail
} // namespace cnf

#endif