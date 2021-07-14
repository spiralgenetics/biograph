
#include <map>
#include <algorithm>
#include <boost/foreach.hpp>
#include <queue>

template<class context>
class astar_state
{
public:
	typedef typename context::dist_t dist_t;
	typedef typename context::location_t location_t;

	astar_state(const context& ctx, const location_t& start, const location_t& goal, const dist_t& max_dist)
		: m_ctx(ctx)
		, m_start(start)
		, m_goal(goal)
		, m_max_dist(max_dist)
	{}

  dist_t run() {
    if (m_start == m_goal) return dist_t();
    m_back.insert(std::make_pair(m_start, m_start));
    proc_location(m_start, dist_t(), 0);
    while (!m_queue.empty()) {
      auto front = m_queue.top();
      if (front.first.dist >= m_max_dist) {
        break;
      }
      m_queue.pop();

      if (m_back.find(front.second.cur) == m_back.end()) {
        m_back.insert(std::make_pair(front.second.cur, front.second.prev));
        if (front.second.cur == m_goal) {
          return front.second.dist;
        }
        proc_location(front.second.cur, front.second.dist,
                      front.first.generation);
      }
    }
    return m_max_dist;
  }

	void get_path(std::vector<location_t>& out)
	{
		location_t cur = m_goal;
		while(cur != m_start)
		{
			if (m_back.find(cur) == m_back.end())
				break;
			out.push_back(cur);
			cur = m_back.find(cur)->second;
		}
		out.push_back(m_start);
		std::reverse(out.begin(), out.end());
	}

private:
	void proc_location(const location_t& loc, const dist_t& dist, size_t generation)
	{
		typedef std::pair<dist_t, location_t> pair_t;
		BOOST_FOREACH(const pair_t& p, m_ctx.nearby(loc))
		{
			if (m_back.find(p.second) != m_back.end()) continue;
			key_t k;
			k.dist = dist + p.first + m_ctx.estimate(p.second, m_goal);
			if (k.dist >= m_max_dist) continue;
			k.generation = generation + 1;
			value_t v(p.second, loc, dist + p.first);
			m_queue.emplace(k, v);
		}
	}
	
	struct key_t 
	{
		dist_t dist;
		size_t generation;
		bool operator<(const key_t& rhs) const
		{
			if (dist != rhs.dist) return dist < rhs.dist;
			return generation < rhs.generation;
		}
	};

	struct value_t
	{
		value_t(const location_t& _cur, const location_t& _prev, const dist_t& _dist)
			: cur(_cur), prev(_prev), dist(_dist) {}	
		location_t cur;
		location_t prev;
		dist_t dist;
	};

	const context& m_ctx;
	const location_t& m_start;
	const location_t& m_goal;
	const dist_t& m_max_dist;
	
  struct key_less {
    bool operator()(const std::pair<key_t, value_t>& p1,
                    const std::pair<key_t, value_t>& p2) const {
      return p2.first < p1.first;
    }
  };

	typedef std::priority_queue<std::pair<key_t, value_t>,
                                std::vector<std::pair<key_t, value_t>>, key_less>
        queue_t;
	queue_t m_queue;
	std::map<location_t, location_t> m_back;
};

