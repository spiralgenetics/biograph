#include "modules/variants/reversable_tracer.h"
#include "modules/variants/align.h"
#include "modules/variants/anchor_drop.h"
#include "modules/variants/pair_counter.h"
#include "modules/variants/sort.h"

#include "modules/bio_base/readmap.h"

namespace variants {

reversable_tracer::reversable_tracer(bool rev_comp, const assemble_options& options)
    : m_rev_comp(rev_comp), m_options(options) {
  CHECK(m_options.scaffold);
  CHECK(!m_options.scaffold->empty());
  m_ref_end_pos = m_options.scaffold->end_pos();

  if (m_rev_comp) {
    m_rev_scaffold = m_options.scaffold->rev_comp();
    m_options.scaffold = &m_rev_scaffold;
    if (m_options.report_half_aligned_func) {
      auto orig_report_half_aligned = m_options.report_half_aligned_func;
      m_options.report_half_aligned_func = [this, orig_report_half_aligned](
          const half_aligned_assembly& ha) {
        CHECK(!ha.right_anchor) << "Forward-only tracer reported a reversed half aligned assembly?";
        orig_report_half_aligned(reverse_half_aligned(ha, m_options.readmap, m_ref_end_pos));
      };
    }
    if (m_options.report_anchor_drop_func) {
      auto orig_report_anchor_drop = m_options.report_anchor_drop_func;
      m_options.report_anchor_drop_func = [this, orig_report_anchor_drop](
          const assembly& a, bool right_anchor) {
        assembly rc_a = a;
        reverse_assembly_in_place(&rc_a, m_options.readmap, m_ref_end_pos);
        orig_report_anchor_drop(rc_a, !right_anchor);
      };
    }
    CHECK(!m_rev_scaffold.empty());
  }

  if (m_options.use_bidir_tracer) {
    m_bidir_tracer.emplace(m_options);
  } else if (m_options.use_pop_tracer) {
    m_pop_tracer.emplace(m_options);
  } else {
    m_tracer.emplace(m_options);
  }
}

std::function<void(const assembly&, bool)>
reversable_tracer::wrap_report_anchor_drop_for_pop_tracer(
    const std::function<void(const assembly&, bool)>& orig_f) {
  if (!m_pop_tracer) {
    return orig_f;
  }

  if (m_rev_comp) {
    return [this, orig_f](const assembly& a, bool right_anchor) {
      assembly rc_a = a;
      reverse_assembly_in_place(&rc_a, m_options.readmap, m_ref_end_pos);
      m_pop_tracer->add_anchor_drop(rc_a, !right_anchor);
      if (orig_f) {
        orig_f(a, right_anchor);
      }
    };
  } else {
    return [this, orig_f](const assembly& a, bool right_anchor) {
      m_pop_tracer->add_anchor_drop(a, right_anchor);
      if (orig_f) {
        orig_f(a, right_anchor);
      }
    };
  }
}

assemble_stats reversable_tracer::assemble(assemble_pipeline_interface* output,
                                           progress_handler_t progress) {
  return assemble(0, m_ref_end_pos, output, progress);
}

assemble_stats reversable_tracer::assemble(aoffset_t start_offset, aoffset_t limit_offset,
                                           assemble_pipeline_interface* output,
                                           progress_handler_t progress) {
  if (limit_offset > m_ref_end_pos) {
    limit_offset = m_ref_end_pos;
  }
  pipeline_step_t tracer_output;
  if (m_rev_comp) {
    if (m_options.only_trace_forward) {
      return assemble_stats();
    }
    tracer_output = make_unique<assemble_lambda_output>(
        [this, output](assembly_ptr a) {
          reverse_assembly_in_place(a.get(), m_options.readmap, m_ref_end_pos);
          output->add(std::move(a));
        },
        "reversing_tracer_rev_output");
  } else {
    tracer_output = make_unique<assemble_lambda_output>(
        [output](assembly_ptr a) {
          output->add(std::move(a));
        },
        "reversing_tracer_fwd_output");
  }

  if (m_options.use_bidir_tracer) {
    if (m_rev_comp) {
      return assemble_stats();
    }
    auto s = make_unique<sorter>(assembly::left_offset_less_than, std::move(tracer_output));
    m_bidir_tracer->add_reference(start_offset, limit_offset);
    m_bidir_tracer->assemble(s.get(), progress);
    return assemble_stats();
  } else if (m_options.use_pop_tracer) {
    auto pc = make_unique<pair_counter>(m_options, std::move(tracer_output));
    auto s = make_unique<sorter>(assembly::left_offset_less_than, std::move(pc));
    if (m_rev_comp) {
      m_pop_tracer->add_reference(m_ref_end_pos - limit_offset, m_ref_end_pos - start_offset);
    } else {
      m_pop_tracer->add_reference(start_offset, limit_offset);
    }
    m_pop_tracer->assemble(s.get());
    return assemble_stats();
  } else {
    auto pc = make_unique<pair_counter>(m_options, std::move(tracer_output));
    auto align = make_unique<anchor_dropper>(m_options, std::move(pc));

    if (m_rev_comp) {
      return m_tracer->assemble(m_ref_end_pos - limit_offset, m_ref_end_pos - start_offset,
                                align.get(), progress);
    } else {
      return m_tracer->assemble(start_offset, limit_offset, align.get(), progress);
    }
  }
}

void reversable_tracer::add_approx_read(uint32_t read_id, aoffset_t start_limit, aoffset_t end_limit,
                                 bool rev_comp) {
  if (!m_pop_tracer) {
    return;
  }
  if (m_rev_comp != rev_comp) {
    read_id = m_options.readmap->get_rev_comp(read_id);
    std::swap(start_limit, end_limit);
    start_limit = m_ref_end_pos - start_limit;
    end_limit = m_ref_end_pos - end_limit;
  }

  m_pop_tracer->add_read(read_id, start_limit, end_limit);
}

}  // namespace variants
