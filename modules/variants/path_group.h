#pragma once

#include "modules/bio_base/seqset.h"

#include <unordered_set>
#include <vector>

namespace variants {

// Path grouper allows tracing through a number of paths on a seqset at once.
class path_group {
 public:
  class listener {
   public:
    virtual void on_seqset_entry(const seqset_range& r, path_group* pg) = 0;
    // virtual int distance_limit_for_range(const seqset_range& r) = 0;

    virtual void on_path_trim(size_t paths) {}

   protected:
    listener() = default;
    ~listener() = default;
  };
  class distant_object {
   public:
    virtual ~distant_object() {}

   protected:
    distant_object() = default;
  };
  class dobj_visitor {
   public:
    virtual void visit(distant_object*, int distance) = 0;

   protected:
    dobj_visitor() = default;
    ~dobj_visitor() = default;
  };

  // Minimum overlap specifies the minimum amount of overlaps on reads
  // if we have more than one path.  Higher values will give better
  // performance because of low branching; lower values will search
  // for paths better.
  path_group(const seqset_range& pos, unsigned min_overlap, listener* l);
  path_group(const path_group&) = delete;
  virtual ~path_group();

  void add_distant_object(std::shared_ptr<distant_object> dobj, int size);
  void visit_distant_objects(const seqset_range& r, dobj_visitor& v);

  void add_sequence(const dna_slice& seq);
  void add_base(dna_base b);

  // Joins two path groups.
  void join(std::unique_ptr<path_group> rhs);

  std::unique_ptr<path_group> split();
  bool empty() const { return m_cur.empty(); }
  size_t size() const { return m_cur.size(); }
  void flush();

  void dump_debug_state();
  // Maximum number of coverage paths to track in parallel.  If zero, unlimited.
  void set_max_size(size_t new_max_size) { m_max_size = new_max_size; }

 protected:
  path_group() = default;
  // Joins from another path group; rhs should not be used after this.
  void join_from(path_group* rhs);
  // Splits the current path group to a freshly initialized empty path group.
  void split_into(path_group* result);

  void reset(const seqset_range& r);

 private:
  struct dobj_set {
    using vec_t = std::vector<std::shared_ptr<distant_object>>;

   public:
    dobj_set() = default;
    // Preserve insert order, but deduplicate.  TODO(nils): Can we
    // just turn path_objects_t into a set<std::pair<int,
    // shared_ptr<dobj>>> when we don't care about exactly reproducing
    // old results?
    ~dobj_set();

    // Adds an individual dobj to this dobj set.
    void add(const std::shared_ptr<distant_object>& dobj);

    // Consumes the given dobj_set and adds all its dobjs into this
    // dobj_set.
    void add(dobj_set&& dobjs);

    using iterator = vec_t::iterator;
    using const_iterator = vec_t::const_iterator;
    iterator begin() { return m_dobjs.begin(); }
    iterator end() { return m_dobjs.end(); }
    const_iterator begin() const { return m_dobjs.begin(); }
    const_iterator end() const { return m_dobjs.end(); }
    bool empty() const { return m_dobjs.empty(); }

   private:
    vec_t m_dobjs;
  };
  using path_objects_t = std::map<int, dobj_set>;
  struct path {
    path_objects_t objects;
    // distance to be added to the per-object distance.
    int unprop_distance = 0;
  };
  struct cur_comparer {
    bool operator()(const seqset_range& lhs, const seqset_range& rhs) const {
      if (lhs.size() != rhs.size()) {
        // Sort longest first, so they get discarded last when we're
        // down to a single branch.
        return lhs.size() > rhs.size();
      }
      return lhs.begin() < rhs.begin();
    }
  };

  void add_path(const seqset_range& pos, path head);
  path merge_paths(const seqset_range& r, path path1, path path2);
  void flush_path(const seqset_range& r, path* head);
  void trim_cur();

  static constexpr size_t k_default_max_size = 64;
  size_t m_max_size = k_default_max_size;

  std::map<seqset_range, path, cur_comparer> m_cur;
  unsigned m_cur_size = 0;

  unsigned m_min_overlap = 0;
  listener* m_listener = nullptr;
};

}  // namespace variants
