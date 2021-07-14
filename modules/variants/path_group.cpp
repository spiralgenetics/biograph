#include "modules/variants/path_group.h"

#include <queue>
#include <unordered_set>

namespace variants {

constexpr bool k_pg_dbg = false;

constexpr size_t path_group::k_default_max_size;

void path_group::dump_debug_state() {
  std::cout << "Path group " << (void*)this << " has " << m_cur.size() << " heads:\n";
  for (auto& cur : m_cur) {
    std::cout << " " << cur.first.sequence().rev_comp() << ":";
    for (const auto& dobjs : cur.second.objects) {
      for (const auto& dobj : dobjs.second) {
        std::cout << " (" << dobj << " + " << cur.second.unprop_distance << ")";
      }
    }
    std::cout << "\n";
  }
}

path_group::path_group(const seqset_range& pos, unsigned min_overlap, listener* l)
    : m_min_overlap(min_overlap), m_listener(l) {
  CHECK_GT(m_max_size, 0);
  add_path(pos, path{});
}

void path_group::reset(const seqset_range& pos) {
  m_cur.clear();
  m_cur_size = 0;
  add_path(pos, path{});
}

void path_group::add_path(const seqset_range& pos, path_group::path p) {
  auto it = m_cur.find(pos);

  if (it == m_cur.end()) {
    CHECK(m_cur.insert(std::make_pair(pos, std::move(p))).second);
    ++m_cur_size;
  } else {
    path merged = merge_paths(pos, std::move(p), std::move(it->second));
    it->second = std::move(merged);
  }
}

void path_group::flush() {
  for (auto& cur : m_cur) {
    flush_path(cur.first, &cur.second);
  }
}

path_group::path path_group::merge_paths(const seqset_range& r, path path1, path path2) {
  if (k_pg_dbg) {
    std::cout << "Merging paths\n";
  }
  if (path1.unprop_distance != path2.unprop_distance) {
    if (path1.unprop_distance) {
      flush_path(r, &path1);
    }
    if (path2.unprop_distance) {
      flush_path(r, &path2);
    }
  }
  CHECK_EQ(path1.unprop_distance, path2.unprop_distance);

  for (auto& obj : path2.objects) {
    path1.objects[obj.first].add(std::move(obj.second));
  }
  return path1;
}

void path_group::flush_path(const seqset_range& r, path* p) {
  if (!p->unprop_distance) {
    return;
  }
  path_objects_t old_objects;
  std::swap(p->objects, old_objects);
  for (auto& obj : old_objects) {
    int adjusted_dist = obj.first + p->unprop_distance;
    if (adjusted_dist <= int(r.size())) {
      p->objects.insert(std::make_pair(adjusted_dist, std::move(obj.second)));
    }
  }
  p->unprop_distance = 0;
}

void path_group::add_distant_object(std::shared_ptr<distant_object> dobj, int size) {
  CHECK(dobj);
  for (auto& cur : m_cur) {
    auto& p = cur.second;
    int adjusted_distance = -(p.unprop_distance + size);
    p.objects[adjusted_distance].add(dobj);
  }
}

void path_group::add_sequence(const dna_slice& seq) {
  //  std::cout << "Adding to path group: " << seq << "\n";
  for (dna_base b : seq) {
    add_base(b);
  }
}

void path_group::add_base(dna_base b) {
  std::map<seqset_range, path, cur_comparer> old_cur;

  std::swap(old_cur, m_cur);
  m_cur_size = 0;

  for (auto& cur : old_cur) {
    const auto& r = cur.first;
    path& p = cur.second;
    ++p.unprop_distance;

    seqset_range new_r = r.push_front_drop(b.complement(), 0);
    add_path(new_r, std::move(p));
  }

  trim_cur();

  for (auto& cur : m_cur) {
    const auto& r = cur.first;
    auto& p = cur.second;
    CHECK(r.valid());
    if (!p.objects.empty()) {
      int adjusted_dist = p.objects.begin()->first + p.unprop_distance;
      if (adjusted_dist > int(r.size())) {
        if (k_pg_dbg) {
          std::cout << "Flushing adjusted distance " << adjusted_dist << " on path " << r.sequence()
                    << "\n";
        }
        flush_path(r, &p);
      }
    }
    if (k_pg_dbg) {
      std::cout << "path group " << cur.first.sequence() << " is a seqset entry, for pg:\n";
      dump_debug_state();
    }
    m_listener->on_seqset_entry(r, this);
  }
}

void path_group::trim_cur() {
  if (m_max_size == 0) {
    // Unlimited
    return;
  }
  if (m_cur_size <= m_max_size) {
    return;
  }
  m_listener->on_path_trim(m_cur_size);
  while (m_cur_size > m_max_size) {
    auto last = m_cur.end();
    CHECK(last != m_cur.begin());
    --last;
    m_cur.erase(last);
    --m_cur_size;
  }
}

void path_group::visit_distant_objects(const seqset_range& r, dobj_visitor& v) {
  std::set<const distant_object*> notify_done;

  auto it = m_cur.find(r);
  CHECK(it != m_cur.end());
  const auto& p = it->second;
  if (k_pg_dbg) {
    std::cout << "Visiting distant objects " << r.sequence() << " on pg:\n";
    dump_debug_state();
  }

  for (const auto& o : p.objects) {
    int adjusted_dist = o.first + p.unprop_distance;
    if (k_pg_dbg) {
      std::cout << "adjusted dist = " << adjusted_dist << "\n";
    }
    if (adjusted_dist > int(r.size())) {
      break;
    }
    for (const auto& dobj : o.second) {
      if (notify_done.insert(dobj.get()).second) {
        if (k_pg_dbg) {
          std::cout << "Notifying\n";
        }
        v.visit(dobj.get(), adjusted_dist);
      }
    }
  }
}

void path_group::join(std::unique_ptr<path_group> rhs) { join_from(rhs.get()); }

void path_group::join_from(path_group* rhs) {
  //  std::cout << "Joining paths, lhs size = " << lhs->m_cur.size()
  //            << " rhs size = " << rhs->m_cur.size() << "\n";
  CHECK_EQ(m_listener, rhs->m_listener);
  CHECK_EQ(m_min_overlap, rhs->m_min_overlap);
  CHECK_EQ(m_cur_size, m_cur.size());
  CHECK_EQ(rhs->m_cur_size, rhs->m_cur.size());
  CHECK_EQ(m_max_size, rhs->m_max_size);
  for (auto& cur : rhs->m_cur) {
    add_path(cur.first, std::move(cur.second));
  }
  CHECK_EQ(m_cur_size, m_cur.size());
  rhs->m_cur.clear();
  rhs->m_cur_size = 0;
}

std::unique_ptr<path_group> path_group::split() {
  std::unique_ptr<path_group> result(new path_group);
  split_into(result.get());
  return result;
}

void path_group::split_into(path_group* result) {
  result->m_min_overlap = m_min_overlap;
  result->m_listener = m_listener;
  result->m_cur = m_cur;
  result->m_cur_size = m_cur_size;
  result->m_max_size = m_max_size;
}

path_group::~path_group() {
  //  std::cout << "Destroying path group with " << m_cur.size() << " paths\n";
  flush();
  CHECK_EQ(m_cur.size(), m_cur_size);
}

void path_group::dobj_set::add(const std::shared_ptr<distant_object>& dobj) {
  for (auto it = m_dobjs.begin(); it != m_dobjs.end(); ++it) {
    if (*it == dobj) {
      // Erase any previous occurrences of this dobj, since we're readding
      // it to the end.
      m_dobjs.erase(it);
      break;
    }
  }
  m_dobjs.push_back(dobj);
}

void path_group::dobj_set::add(dobj_set&& dobjs_to_add) {
  if (dobjs_to_add.empty()) {
    return;
  }
  if (m_dobjs.empty()) {
    m_dobjs = std::move(dobjs_to_add.m_dobjs);
    return;
  }
  std::unordered_set<distant_object*> to_erase;
  for (const auto& dobj : dobjs_to_add.m_dobjs) {
    to_erase.insert(dobj.get());
  }
  for (auto it = m_dobjs.begin(); it != m_dobjs.end();) {
    auto to_erase_it = to_erase.find(it->get());
    if (to_erase_it == to_erase.end()) {
      ++it;
    } else {
      it = m_dobjs.erase(it);
      to_erase.erase(to_erase_it);
    }
  }

  for (auto& dobj : dobjs_to_add.m_dobjs) {
    m_dobjs.push_back(std::move(dobj));
  }
}

path_group::dobj_set::~dobj_set() {
  // Explicit destruction order of elements; see TODO in dobj_set declaration.
  while (!m_dobjs.empty()) {
    m_dobjs.pop_back();
  }
}

}  // namespace variants
