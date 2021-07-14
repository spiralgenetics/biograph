#include "modules/variants/trace_ref.h"

#include "modules/bio_base/readmap.h"
#include "modules/io/parallel.h"
#include "modules/variants/reversable_tracer.h"

namespace variants {

trace_ref::trace_ref(const assemble_options& options, scaffold_pipeline_interface* output_f)
    : m_options(options), m_output_f(output_f) {
  if (m_options.use_pop_tracer || m_options.pop_trace_anchor_drop || m_options.use_bidir_tracer) {
    CHECK(options.readmap);
    options.readmap->calc_read_len_limits_if_needed();
  }

  CHECK(!options.scaffold) << "Scaffold should not already be provided";
  CHECK(options.ref);
}

void trace_ref::add_scaffold(const std::string& scaffold_name) {
  add_scaffold_range(scaffold_name, 0, std::numeric_limits<size_t>::max());
}

void trace_ref::add_scaffold_range(const std::string& scaffold_name, size_t start, size_t limit) {
  CHECK_GE(limit, start);

  bool skip_fwd = m_options.skip_push_trace_fwd;
  bool skip_rev = m_options.skip_push_trace_rev;
  CHECK(!(skip_fwd && skip_rev)) << "No tracing specified?";

  std::shared_ptr<scaffold> s = get_scaffold(scaffold_name);
  CHECK(s);
  assemble_options opts = m_options;
  opts.scaffold = s.get();
  opts.scaffold_name = scaffold_name;

  aoffset_t aostart = start;
  aoffset_t aolimit = s->end_pos();
  if (limit < size_t(aolimit)) {
    aolimit = limit;
  }
  CHECK_GE(aolimit, aostart);

  while (aostart < aolimit) {
    std::unique_ptr<work_info> w = make_unique<work_info>();

    w->p = m_output_f->pipeline_for_scaffold(opts, scaffold_name);
    w->start = aostart;
    w->limit = std::min<aoffset_t>(aolimit, m_options.scaffold_split_size + aostart);
    w->scaffold_name = scaffold_name;
    w->s = s;
    w->skip_fwd = skip_fwd;
    w->skip_rev = skip_rev;

    if (m_options.use_bidir_tracer) {
      aostart = std::min<aoffset_t>(w->limit, aostart + aoffset_t(m_options.scaffold_split_size) -
                                                  m_options.read_ahead_distance);
    } else {
      aostart = w->limit;
    }

    m_work.push_back(std::move(w));
  }
}

assemble_stats trace_ref::assemble(progress_handler_t progress) {
  // std::cout << m_work.size() << " work items\n";

  auto atoi_or_maxint = [](const std::string& s) -> int {
    try {
      return std::stoi(s);
    } catch (const std::invalid_argument&) {
    } catch (const std::out_of_range&) {
    }
    return std::numeric_limits<int>::max();
  };

  std::sort(m_work.begin(), m_work.end(),
            [&](const std::unique_ptr<work_info>& awork, const std::unique_ptr<work_info>& bwork) {
              const work_info& a = *awork;
              const work_info& b = *bwork;
              int an = atoi_or_maxint(a.scaffold_name);
              int bn = atoi_or_maxint(b.scaffold_name);
              if (an != bn) {
                // Numeric scaffolds go in order,
                // non-numeric scaffolds go last.
                return an < bn;
              }
              // Otherwise, shorter scaffold names go first (so e.g. "X" and "Y" go before
              // "hs37d5").
              int as = a.scaffold_name.size();
              int bs = b.scaffold_name.size();
              if (as != bs) {
                return as < bs;
              }
              // Otherwise, group the fwd and reverse versions of chunks together.
              if (a.scaffold_name != b.scaffold_name) {
                return b.scaffold_name < a.scaffold_name;
              }
              return a.start < b.start;
            });

  assemble_stats tot_st;
  std::mutex mu;
  parallel_for(
      0, m_work.size(),
      [&](size_t idx) {
        if (m_aborted) {
          return;
        }
        auto& w = m_work[idx];

        assemble_stats st = execute_work(std::move(w));

        std::lock_guard<std::mutex> l(mu);
        tot_st += st;
      },
      progress);
  if (m_aborted) {
    for (auto& work : m_work) {
      abort_work(std::move(work));
    }
  }
  m_work.clear();
  return tot_st;
}

void trace_ref::abort_work(std::unique_ptr<work_info> w) {
  if (!w) {
    return;
  }
  w->pop.reset();
  w->rc_pop.reset();
  w->p.reset();
  w->s.reset();
  if (g_verbose_trace_work) {
    std::lock_guard<std::mutex> l(g_in_progress_mu);
    SPLOG("ABORT: %s", w->to_string().c_str());
  }
}

assemble_stats trace_ref::execute_work(std::unique_ptr<work_info> w) const {
  CHECK(w->s);
  assemble_stats st;

  if (w->s->empty()) {
    return st;
  }

  assemble_options opts = m_options;
  opts.scaffold_name = w->scaffold_name;
  opts.scaffold = w->s.get();
  CHECK(opts.scaffold);

  assemble_options pop_opts = opts;
  pop_opts.use_pop_tracer = true;
  pipeline_step_t pop_out;

  if (m_options.pop_trace_anchor_drop && !m_options.use_pop_tracer) {
    pop_out = w->p->make_parallel_input();
    auto report_f = m_options.report_anchor_drop_func;
    CHECK(!(m_options.skip_pop_trace_rev && m_options.skip_pop_trace_fwd))
        << "Pop tracing specified, but skipping both directions?";
    if (!m_options.skip_pop_trace_fwd) {
      w->pop = make_unique<reversable_tracer>(false /* not rev comp */, pop_opts);
      report_f = w->pop->wrap_report_anchor_drop_for_pop_tracer(report_f);
    }

    if (!m_options.skip_pop_trace_rev) {
      w->rc_pop = make_unique<reversable_tracer>(true /* rev comp */, pop_opts);
      report_f = w->rc_pop->wrap_report_anchor_drop_for_pop_tracer(report_f);
    }

    w->report_anchor_drop_func = report_f;
  }

  assemble_options push_opts = opts;
  if (w->report_anchor_drop_func) {
    push_opts.report_anchor_drop_func = w->report_anchor_drop_func;
  }

  if (!w->skip_fwd) {
    note_work_start(w, "push-fwd");
    st += execute_work_direction(w.get(), false /* forward */, push_opts);
    note_work_finish(w, "push-fwd");
  }
  if (!w->skip_rev) {
    note_work_start(w, "push-rev");
    st += execute_work_direction(w.get(), true /* reverse */, push_opts);
    note_work_finish(w, "push-rev");
  }

  if (w->pop || w->rc_pop) {
    if (w->pop) {
      note_work_start(w, "pop-fwd");
      st += w->pop->assemble(w->start, w->limit, pop_out.get());
      w->pop.reset();
      note_work_finish(w, "pop-fwd");
    }

    if (w->rc_pop) {
      note_work_start(w, "pop-rev");
      st += w->rc_pop->assemble(w->start, w->limit, pop_out.get());
      w->rc_pop.reset();
      note_work_finish(w, "pop-rev");
    }
    pop_out.reset();
  } else {
    CHECK(!pop_out);
  }

  note_work_start(w, "flush");
  w->p.reset();
  w->s.reset();
  note_work_finish(w, "flush");
  w.reset();

  return st;
}

assemble_stats trace_ref::execute_work_direction(work_info* w, bool rev_comp,
                                                 const assemble_options& opts) const {
  assemble_stats st;
  auto start_time = std::chrono::high_resolution_clock::now();
  reversable_tracer t(rev_comp, opts);
  CHECK(w->p);
  pipeline_step_t out = w->p->make_parallel_input();
  if (w->pop || w->rc_pop) {
    auto new_out = make_unique<assemble_lambda_copy>(
        [&](const assembly& a) {
          if (a.matches_reference) {
            return;
          }
          for (uint32_t rc_read_id : a.rc_read_ids) {
            if (!m_options.readmap->has_mate(rc_read_id)) {
              continue;
            }
            uint32_t read_id = m_options.readmap->get_rev_comp(rc_read_id);
            aoffset_t mate_start_limit = a.left_offset;
            aoffset_t mate_end_limit = a.right_offset;
            if (m_options.readmap->get_is_forward(rc_read_id) ==
                m_options.forward_pairs_face_inward) {
              if (mate_start_limit <= aoffset_t(m_options.max_pair_distance)) {
                mate_start_limit = 0;
              } else {
                mate_start_limit -= m_options.max_pair_distance;
              }
            } else {
              mate_end_limit += m_options.max_pair_distance;
            }
            uint32_t mate_id = m_options.readmap->get_mate(rc_read_id);

            if (w->pop) {
              w->pop->add_approx_read(read_id, mate_start_limit, mate_end_limit, rev_comp);
              w->pop->add_approx_read(mate_id, mate_start_limit, mate_end_limit, rev_comp);
            }
            if (w->rc_pop) {
              w->rc_pop->add_approx_read(read_id, mate_start_limit, mate_end_limit, rev_comp);
              w->rc_pop->add_approx_read(mate_id, mate_start_limit, mate_end_limit, rev_comp);
            }
          }
        },
        std::move(out), "save_reads_for_pop");
    out = std::move(new_out);
  }
  st = t.assemble(w->start, w->limit, out.get());
  if (m_options.report_chunk_stats_func) {
    auto end_time = std::chrono::high_resolution_clock::now();
    double duration = std::chrono::duration<double>(end_time - start_time).count();
    m_options.report_chunk_stats_func(w->scaffold_name, w->start, w->limit, rev_comp, duration, st);
  }
  return st;
}

void trace_ref::add_entire_reference() {
  for (const auto& scaffold_name : m_options.ref->get_assembly().scaffold_order) {
    add_scaffold(scaffold_name);
  }
}

std::shared_ptr<scaffold> trace_ref::get_scaffold(const std::string& scaffold_name) const {
  // Reuse existing if possible.
  for (const auto& w : m_work) {
    if (w->scaffold_name == scaffold_name) {
      return w->s;
    }
  }

  return std::make_shared<scaffold>(ref_to_scaffold(m_options.ref, scaffold_name));
}

scaffold trace_ref::ref_to_scaffold(const reference* ref, const std::string& scaffold_name) {
  scaffold result;

  for (const auto& sc : ref->get_assembly().supercontigs) {
    if (sc.scaffold_name == scaffold_name) {
      result.add(sc.offset,
                 dna_slice(ref->get_dna(sc.tot_offset), ref->get_dna(sc.tot_offset + sc.len)));
    }
  }
  return result;
}

std::mutex trace_ref::g_in_progress_mu;
std::map<trace_ref::in_progress_key_t, time_t /* start time */> trace_ref::g_in_progress;
bool trace_ref::g_verbose_trace_work = false;

void trace_ref::note_work_start(const std::unique_ptr<work_info>& w, const std::string& work_desc) {
  bool did_insert;
  {
    std::lock_guard<std::mutex> l(g_in_progress_mu);
    in_progress_key_t k = {w.get(), work_desc};
    did_insert = g_in_progress.emplace(k, time(0)).second;
  }
  if (!did_insert || g_verbose_trace_work) {
    SPLOG("START: %s %s", w->to_string().c_str(), work_desc.c_str());
    SPLOG("TRACES: %s", work_in_progress().c_str());
  }
  CHECK(did_insert) << "Duplicate work note start? work:" << w->to_string()
                    << " desc: " << work_desc;
}

void trace_ref::note_work_finish(const std::unique_ptr<work_info>& w,
                                 const std::string& work_desc) {
  bool did_find;
  time_t start_time = 0;
  {
    std::lock_guard<std::mutex> l(g_in_progress_mu);
    in_progress_key_t k = {w.get(), work_desc};
    auto it = g_in_progress.find(k);
    did_find = it != g_in_progress.end();
    if (did_find) {
      start_time = it->second;
      g_in_progress.erase(it);
    }
  }

  if (g_verbose_trace_work || !did_find) {
    time_t now = time(0);
    SPLOG("FINISH: %s %s (%lds)", w->to_string().c_str(), work_desc.c_str(), now - start_time);
    SPLOG("TRACES: %s", work_in_progress().c_str());
  }
  CHECK(did_find) << "Note work finish for missing work? work: " << w->to_string()
                  << " desc: " << work_desc << " in progress: " << work_in_progress();
}

std::string trace_ref::work_in_progress() {
  std::stringstream result;
  std::lock_guard<std::mutex> l(g_in_progress_mu);
  time_t now = time(0);
  result << g_in_progress.size() << " in progress:";
  for (const auto& ip : g_in_progress) {
    const auto& w = ip.first.first;
    const auto& work_desc = ip.first.second;
    time_t start = ip.second;
    int duration = now - start;

    result << " " << w->to_string() << " " << work_desc << "(" << duration << "s)";
  }
  return result.str();
}

trace_ref::work_info::~work_info() { CHECK(!p) << "Pipeline should have been flushed."; }

bool trace_ref::empty() const { return m_work.empty(); }

void trace_ref::abort_trace() { m_aborted = true; }

trace_ref::~trace_ref() {
  if (m_aborted) {
    for (auto& work : m_work) {
      if (work) {
        abort_work(std::move(work));
      }
    }
  }
  m_work.clear();
}

}  // namespace variants
