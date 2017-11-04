#include "cnf/cnf.hpp"
#include "crypto/sha256.hxx"
#include "cnf/multivar.hpp"

#include <iostream>
#include <cassert>
#include <cstring>
#include <vector>
#include <string>
#include <ctime>

using namespace std;
using namespace crypto;

inline uint32_t bswap(uint32_t x)
{
    x = (x & 0xff00ff) << 8 | ((x >> 8) & 0xff00ff);
    return (x >> 16) | (x << 16);
}

vector<uint32_t> sha256(const char *str) {
    vector<uint32_t> h(8);
    sha256(h.data(), str, strlen(str));
    return h;
}

vector<uint32_t> sha256(const std::string &str) {
    vector<uint32_t> h(8);
    sha256(h.data(), str.data(), str.size());
    return h;
}

string hex(const vector<uint32_t> &data) {
    string res(data.size() * 8, ' ');
    for (size_t i = 0; i < data.size(); ++i) {
        uint32_t x = data[i];
        for (int j = 8; j --> 0; ) {
            int l = x % 16;
            x /= 16;
            res[i * 8 + j] = (l < 10 ? '0' : 'a' - 10) + l;
        }
    }
    return res;
}

void print(cnf::var_pool *pool, const cnf::tree_node &node)
{
	if (node.type() == cnf::tree_node::LEAF) {
		assert(node.key());
		cout << node.key();
	} else {
		std::string sep = node.type() == cnf::tree_node::AND ? " & " : " | ";
		cout << '(';
		bool first = true;
		for (int x : node) {
			if (first)
				first = false;
			else
				cout << sep;
			if (x < 0)
				cout << '~';
			print(pool, pool->node(x));
		}
		cout << ')';
	}
}

void print(const cnf::var &var)
{
	if (var.is_constant())
		cout << (var == true);
	else {
		if (var.id() < 0)
			cout << '~';
		print(var.pool(), var.node());
	}
}

void println(const cnf::var &var)
{
	print(var);
	cout << endl;
}

template <int n>
void print(const cnf::multivar<n> &var)
{
	for (int i = n; i --> 0; ) {
		print(var[i]);
		if (i) cout << ' ';
	}
}

void test_mv()
{
	auto x = cnf::int2multivar(16777215);
	auto y = cnf::int2multivar(257);
	print(x); cout << endl;
	print(y); cout << endl;
	print(x + y); cout << endl;
	print(x ^ y); cout << endl;
}

int main()
{
	cnf::var_pool vp;
	auto x = vp.make("x");
	auto y = vp.make("y");
	auto z = vp.make("z");
	auto c0 = cnf::var(false);
	auto c1 = cnf::var(true);
	println(~x);
	println(y);
	println(z);
	println(~x & ~(~y | z) | x);
	println(x & y & z & ~x);
	println(x ^ y ^ z);
	println(c0);
	println(c1);
	test_mv();
	return 0;
}
