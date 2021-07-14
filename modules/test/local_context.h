#include "modules/mapred/path.h"
#include "modules/mapred/manifest.h"
#include "modules/mapred/task_mgr.h"
#include "modules/mapred/map_reduce_task.h"
#include "modules/mapred/sort_task.h"
#include "modules/mapred/map_task.h"
#include "modules/mapred/taskdb.h"
#include "modules/bio_format/exporter.h"
#include "modules/io/config.h"
#include "modules/io/make_unique.h"

class local_context
{
public:
	local_context(int num_partitions, int goal_size, const path& out_path)
		: m_num_partitions(num_partitions)
		, m_goal_size(goal_size)
		, m_out_path(out_path)
	{}

	manifest map_only(const std::string& map, const std::string& map_param, const manifest& input, bool is_pipe = false)
	{
		std::unique_ptr<map_task> mt = make_unique<map_task>();

		mt->input = input;
		mt->map = map;
		mt->map_param = map_param;
		mt->is_pipe = is_pipe;

		mt->output_goal_size = m_goal_size;

		manifest out;
		m_task_mgr.run_task(out, m_out_path, std::move(mt));
		return out;
	}

	manifest map_reduce(
		const std::string& map, const std::string& map_param,
		const std::string& sort,
		const std::string& reduce, const std::string& reduce_param,
		const manifest& input,
		bool summarize = false)
	{
		std::unique_ptr<map_reduce_task> mr = make_unique<map_reduce_task>();
		mr->input = input;

		mr->map = map;
		mr->map_param = map_param;
		mr->sort = sort;
		mr->reduce = reduce;
		mr->reduce_param = reduce_param;

		mr->is_summary = summarize;
		mr->num_partitions = m_num_partitions;
		mr->temp_goal_size = m_goal_size;
		mr->output_goal_size = m_goal_size;

		manifest out;
		m_task_mgr.run_task(out, m_out_path, std::move(mr));
		return out;
	}

	manifest map_sort_reduce(
		const std::string& map, const std::string& map_param,
		const std::string& in_sort, const std::string& out_sort,
		const std::string& reduce, const std::string& reduce_param,
		const manifest& input,
		bool summarize = false,
		bool is_pipe = false)
	{
		std::unique_ptr<map_task> mt = make_unique<map_task>();
		mt->input = input;
		mt->map = map;
		mt->map_param = map_param;
		mt->sort = in_sort;
		mt->output_goal_size = m_goal_size;
		mt->is_pipe = is_pipe;

		manifest mapped;
		m_task_mgr.run_task(mapped, m_out_path, std::move(mt));
		std::unique_ptr<sort_task> st = make_unique<sort_task>();
		st->input = mapped;
		st->goal_size = m_goal_size;
		manifest sorted;
		m_task_mgr.run_task(sorted, m_out_path, std::move(st));
		SPLOG("Sorted manifest: %s\n", json_serialize(sorted).c_str());
		std::unique_ptr<sorted_reduce_task> sr = make_unique<sorted_reduce_task>();
		sr->input = sorted;
		sr->reduce = reduce;
		sr->reduce_param = reduce_param;
		sr->out_sort = out_sort;
		sr->prereduce_goal_size = m_goal_size;
		sr->goal_size = m_goal_size;
		manifest final;
		m_task_mgr.run_task(final, m_out_path, std::move(sr));
		SPLOG("Reduced manifest: %s\n", json_serialize(final).c_str());
		return final;
	}

private:
	int m_num_partitions;
	int m_goal_size;
	path m_out_path;
	task_mgr_local m_task_mgr;
};

template<class exporter_type>
void simple_export(const std::string& outfile, const manifest& data)
{
	path output_path(outfile);
	std::unique_ptr<writable> output_writable(output_path.write());
	manifest_reader the_manifest_reader(data);
	exporter_type impl(*output_writable);
	impl.export_from(the_manifest_reader);
}
