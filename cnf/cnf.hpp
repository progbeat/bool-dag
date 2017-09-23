#ifndef CNF_CNF_HPP_
#define CNF_CNF_HPP_

#include "detail/cnf.hpp"

#include <unordered_map>
#include <cstring>
#include <cassert>
#include <memory>
#include <vector>
#include <queue>

namespace cnf {

class var_pool;

// (-2, -1, 0, 1) -> (3, 1, 0, 2)
inline int neg2odd(int x) {
	const int inv = x >> (sizeof(int) * 8 - 1);
	return (x ^ inv) * 2 - inv;
}

class tree_node
{
public:
	enum {
		LEAF = 0,
		OR = '|',
		AND = '&',
	};
	int type() const { return type_; }
	const char *key() const { return key_; }
	const int *begin() const { return begin_; }
	const int *end() const { return end_; }
	tree_node()
		: begin_(nullptr)
		, end_(nullptr)
		, key_(nullptr)
		, type_(LEAF)
	{}
private:
	const int *begin_;
	const int *end_;
	const char *key_;
	int type_;
	
	tree_node(const int *data)
	{
		type_ = data[0] & 0xff;
		if (type_ != LEAF) {
			begin_ = data + 1;
			end_ = data + (data[0] >> 8);
			key_ = nullptr;
		} else {
			begin_ = end_ = nullptr;
			key_ = (const char*)(data + 1);
		}
	}
	
	friend class var_pool;
};

class var
{
public:	
	~var();
	
	var operator ~ () const {
		return is_constant()
			? var(!x_)
			: var(~x_, pool_);
	}
	friend var operator & (const var &lhs, const var &rhs);
	friend var operator | (const var &lhs, const var &rhs);
	
	tree_node node() const;
	
	int id() const { return x_; }
	
	bool operator == (bool value) const
	{
		return pool_ == nullptr && (x_ != 0) == value;
	}
	
	var(bool value = false) : x_(value), pool_(nullptr) {}

	var(const var &rhs) { *this = rhs; }
	
	var(var &&rhs) : x_(rhs.x_), pool_(rhs.pool_)
	{
		rhs.pool_ = nullptr;
	}
	
	var& operator = (const var &rhs);
	
	var& operator = (var &&rhs)
	{
		std::swap(x_, rhs.x_);
		std::swap(pool_, rhs.pool_);
		return *this;
	}
	
	bool is_constant() const { return pool_ == nullptr; }
	
	var_pool *pool() const { return pool_; }
private:
	int x_;
	var_pool *pool_;

	var(int, var_pool*);
	
	friend class var_pool;
};

class var_pool
{
public:
	var make(const char *key)
	{
		auto it = keys_map_.find(key);
		return var(
			it == keys_map_.end() ? make_leaf(key) : it->second,
			this
		);
	}
	
	int size() { return static_cast<int>(pool_.size()); }
	
	tree_node node(int i) const { return tree_node(data(i)); }
	
	static var_pool* from_vars() { return nullptr; }
	
	template <class... Args>
	static var_pool* from_vars(const var &v, Args&&... args)
	{
		var_pool* ptr = from_vars(std::forward<Args>(args)...);
		if (v.pool_ == nullptr)
			return ptr;
		if (ptr == nullptr)
			return v.pool_;
		assert(ptr == v.pool_);
		return ptr;
	}
private:
	typedef std::unique_ptr<int[]> data_ptr;
	std::vector<data_ptr> pool_;
	std::priority_queue<int, std::vector<int>, std::greater<int>> holes_;
	std::unordered_map<std::string, int> keys_map_;
	
	static int index(int i) { return i ^ (i >> (sizeof(int) * 8 - 1)); }
	
	int make_empty(int size)
	{
		assert(size > 0);
		int *ptr = new int[size + 1];
		ptr[0] = 0;
		if (holes_.empty()) {
			pool_.emplace_back(ptr);
			return pool_.size() - 1;
		}
		int x = holes_.top();
		holes_.pop();
		pool_[x].reset(ptr);
		return x;
	}
	
	int make_leaf(const char *key)
	{
		auto len = strlen(key);
		int id = make_empty((len + sizeof(int) * 2 - 1) / 4);
		int *ptr = data(id);
		ptr[-1] = 1;
		ptr[0] = tree_node::LEAF;
		memcpy(ptr + 1, key, len + 1);
		return id;
	}
	
	template <class Ptr> int make_inner(int action, Ptr childs, int count)
	{
		assert(count > 0);
		assert(
			action == tree_node::OR ||
			action == tree_node::AND
		);
		int id = make_empty(++count);
		int *ptr = data(id);
		ptr[0] = action | (count << 8);
		for (int i = 1; i < count; ++i) {
			ptr[i] = *childs;
			++childs;
		}
		return id;
	}
	
	int* data(int i) const { return pool_[index(i)].get() + 1; }
	
	void erase(int i)
	{
		i = index(i);
		if (pool_[i] != nullptr) {
			if (i + 1 < pool_.size()) {
				pool_[i].reset();
				holes_.push(i);
			} else {
				pool_.pop_back();
			}
		}
	}
	
	void take(int i, int cnt)
	{
		i = index(i);
		int *d = data(i);
		if ((d[-1] += cnt) == 0) {
			
		}
	}
	
	friend class var;
	friend var operator & (const var &lhs, const var &rhs);
	friend var operator | (const var &lhs, const var &rhs);
};

var operator & (const var &lhs, const var &rhs)
{
	if (lhs == false || rhs == false)
		return var(false);
	std::unordered_set<int> childs;
	auto process = [&childs](const var &v) {
		if (v == true)
			return;
		auto node = v.node();
		if (node.type() == tree_node::LEAF
			|| (v.id() > 0) != (node.type() == tree_node::AND))
		{
			childs.insert(v.id());
		} else {
			int inv = v.id() >> (sizeof(int) * 8 - 1);
			for (int x : node) 
				childs.insert(x ^ inv);
		}
	};
	process(lhs);
	process(rhs);
	if (childs.empty())
		return var(true);
	auto pool = var_pool::from_vars(lhs, rhs);
	if (childs.size() == 1)
		return var(*childs.begin(), pool);
	else if (detail::has_conflict(childs))
		return var(false);
	return var(
		pool->make_inner(tree_node::AND, childs.begin(), childs.size()),
		pool
	);
}

var operator | (const var &lhs, const var &rhs)
{
	if (lhs == true || rhs == true)
		return var(true);
	std::unordered_set<int> childs;
	auto process = [&childs](const var &v) {
		if (v == false)
			return;
		auto node = v.node();
		if (node.type() == tree_node::LEAF
			|| (v.id() > 0) != (node.type() == tree_node::OR))
		{
			childs.insert(v.id());
		} else {
			int inv = v.id() >> (sizeof(int) * 8 - 1);
			for (int x : node) 
				childs.insert(x ^ inv);
		}
	};
	process(lhs);
	process(rhs);
	if (childs.empty())
		return var(false);	
	auto pool = var_pool::from_vars(lhs, rhs);
	if (childs.size() == 1)
		return var(*childs.begin(), pool);
	else if (detail::has_conflict(childs))
		return var(true);
	return var(
		pool->make_inner(tree_node::OR, childs.begin(), childs.size()),
		pool
	);
}

var operator ^ (const var &lhs, const var &rhs)
{
	return ((~lhs) & rhs) | (lhs & (~rhs));
}

inline tree_node var::node() const
{
	return pool_ != nullptr ? pool_->node(x_) : tree_node();
}

var::var(int x, var_pool *pool)
	: x_(x)
	, pool_(pool)
{
	assert(pool_ != nullptr);
	pool_->take(x, 1);
}

var::~var()
{
	if (pool_ != nullptr)
		pool_->take(x_, -1);
}

var& var::operator = (const var &rhs)
{
	x_ = rhs.x_;
	pool_ = rhs.pool_;
	if (pool_)
		pool_->take(x_, 1);
	return *this;
}

} // namespace cnf

#endif