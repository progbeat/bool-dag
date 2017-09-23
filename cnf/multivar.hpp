#ifndef CNF_MULTIVAR_HPP_
#define CNF_MULTIVAR_HPP_

#include "cnf.hpp"

namespace cnf {

template <int n> class multivar
{
public:
	int size() const { return n; }
	
	var& operator[] (int i)
	{
		assert(0 <= i && i < n);
		return vars_[i];
	}
	const var& operator[] (int i) const
	{
		assert(0 <= i && i < n);
		return vars_[i];
	}
	multivar<n> operator << (int sh) const
	{
		if (sh == 0)
			return *this;
		assert(0 < sh && sh < n);
		multivar<n> res;
		for (int i = sh; i < n; ++i)
			res[i] = vars_[i - sh];
		return res;
	}
	multivar<n> operator >> (int sh) const
	{
		if (sh == 0)
			return *this;
		assert(0 < sh && sh < n);
		multivar<n> res;
		for (int i = sh; i < n; ++i)
			res[i - sh] = vars_[i];
		return res;
	}
private:
	var vars_[n];
};

template <int n>
multivar<n> operator | (const multivar<n> &lhs, const multivar<n> &rhs)
{
	multivar<n> r;
	for (int i = 0; i < n; ++i)
		r[i] = lhs[i] | rhs[i];
	return r;
}

template <int n>
multivar<n> operator & (const multivar<n> &lhs, const multivar<n> &rhs)
{
	multivar<n> r;
	for (int i = 0; i < n; ++i)
		r[i] = lhs[i] & rhs[i];
	return r;
}

template <int n>
multivar<n> operator ^ (const multivar<n> &lhs, const multivar<n> &rhs)
{
	multivar<n> r;
	for (int i = 0; i < n; ++i)
		r[i] = lhs[i] ^ rhs[i];
	return r;
}

template <int n>
multivar<n> operator + (const multivar<n> &lhs, const multivar<n> &rhs)
{
	multivar<n> r;
	r[0] = lhs[0] ^ rhs[0];
	var carry;
	for (int i = 0; i < n - 1; ) {
		carry = (lhs[i] & rhs[i]) | (carry & (lhs[i] | rhs[i]));
		++i;
		r[i] = lhs[i] ^ rhs[i] ^ carry;
	}
	return r;
}

template <class X> multivar<sizeof(X) * 8> int2multivar(X x)
{
	static_assert(std::is_integral<X>::value, "X isn't an integral type");
	multivar<sizeof(X) * 8> res;
	for (int i = 0; i < res.size(); ++i)
		res[i] = var(bool((x >> i) & 1));
	return res;
}

} // namespace cnf

#endif