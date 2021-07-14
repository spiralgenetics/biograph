
#include "modules/mapred/mapred_wrapper_task.h"
#include "modules/mapred/run_wrapper_task.h"

REGISTER_TASK(mapred_wrapper_task);

void mapred_wrapper_task::run() {
	if (m_state == 0) {
		size_t goal_size = std::min(size_t(64*1024*1024), input.get_size() / parts);
		std::vector<input_stream_params> params;
		input.split_by_goal_size(params, goal_size);
		for (const auto& ip : params) {
			auto t = make_unique<run_wrapper_task>();
			for (const auto& m : aux_map_inputs) {
				t->inputs.push_back(m);
			}
			manifest m;
			for (const auto& fi : ip.inputs) {
				m.add(fi, 0);
			}
			t->inputs.push_back(m);
			t->params = map_params;
			t->num_outputs = 1;
			m_map_task.push_back(add_subtask(std::move(t)));
		}
		m_state = 1;
	} else if (m_state == 1) {
		auto t = make_unique<run_wrapper_task>();
		for (const auto& m : aux_map_inputs) {
			t->inputs.push_back(m);
		}
		t->params = reduce_params;
		for (const auto& id : m_map_task) {
			std::vector<manifest> out;
			get_output(out, id);
			t->inputs.push_back(out[0]);
		}
		t->params.args.push_back("--in-count");
		t->params.args.push_back(std::to_string(m_map_task.size()));
		t->num_outputs = 1;
		m_reduce_task = add_subtask(std::move(t));
		m_state = 2;
	} else {
		std::vector<manifest> out;
		get_output(out, m_reduce_task);
		set_output(out[0]);
	}
}

