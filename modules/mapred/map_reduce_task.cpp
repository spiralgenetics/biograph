#include "modules/mapred/map_reduce_task.h"
#include "modules/mapred/sort_task.h"
#include "modules/io/make_unique.h"

#include <memory>

REGISTER_TASK(map_reduce_task);

template<class Task>
void set_params(Task* t, map_reduce_task* mr, const manifest& input)
{
	t->input = input;
	if (mr->add_before_reduce.get_num_records() != 0) {
		t->input.add(mr->add_before_reduce);
	}
	t->reduce = mr->reduce;
	t->reduce_param = mr->reduce_param;
	t->is_summary = mr->is_summary;
	t->goal_size = mr->output_goal_size;
	t->update_freq = mr->reduce_update_freq;
}

void map_reduce_task::run()
{
	if (m_state == 0) {
		split_progress(.005, .5);
		std::unique_ptr<map_task> mt = make_unique<map_task>();
		mt->input = input;
		mt->map = map;
		mt->map_param = map_param;
		mt->input_goal_size = input_goal_size;
		mt->output_goal_size = temp_goal_size;
		mt->update_freq = map_update_freq;
		mt->num_partitions = num_partitions;
		mt->sort = sort;
		if (is_summary) {
			mt->reduce = reduce;
			mt->reduce_param = reduce_param;
		}
		m_map_task = add_subtask(std::move(mt));
		m_state = 1;
	}
	else if (m_state == 1) {
		split_progress(0.005, .01);
		manifest map_result;
		get_output(map_result, m_map_task);
		if (is_summary && use_sort)
		{
			std::unique_ptr<sort_task> st = make_unique<sort_task>();
			set_params(st.get(), this, map_result);
			m_reduce_task = add_subtask(std::move(st));
		}
		else if (use_sort && reduce == "identity")
		{
			std::unique_ptr<sort_task> st = make_unique<sort_task>();
			st->input = map_result;
			m_reduce_task = add_subtask(std::move(st));
		}
		else
		{
			std::unique_ptr<reduce_task> rt = make_unique<reduce_task>();
			set_params(rt.get(), this, map_result);
			rt->mp_goal_size = mp_goal_size;
			rt->goal_size = output_goal_size;
			rt->post_sort = post_sort;
			m_reduce_task = add_subtask(std::move(rt));
		}
		m_state = 2;
	}
	else {
		manifest reduce_result;
		get_output(reduce_result, m_reduce_task);
		set_output(reduce_result);
	}
}
