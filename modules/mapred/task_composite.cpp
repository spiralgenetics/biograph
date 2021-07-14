#include "modules/mapred/task_composite.h"
#include "modules/mapred/task_tree.h"
#include "modules/io/log.h"

bool starts_with(const std::string& prefix, const std::string& source)
{
	return source.find(prefix) == 0UL;
}

void for_each_child_in(
		const std::vector<std::string>& subtasks,
		component_handler_t on_composite,
		component_handler_t on_composite_subtask,
		component_handler_t on_leaf )
{
	for ( auto it = subtasks.cbegin(); it != subtasks.cend(); )
	{
		const auto& a_task = *it;
		if( starts_with("parallel", a_task) || starts_with("serial", a_task) )
		{
			on_composite( it );
			while( (++it != subtasks.cend()) && (*it != a_task) )
			{
				on_composite_subtask(it);
			}


			// the following 2 throws are assertions.
			// If they happen, that means something is logically inconsistent in this class
			// regardless of its input.
			if ( it == subtasks.cend() )
			{
				throw io_exception("invalid subtask group");
			}
			if ( *it != a_task )
			{
				throw io_exception("invalid subtask group: expected " + a_task + " got: " + (*it));
			}
		}
		else
		{
			on_leaf( it );
		}
		it++;
	}
}

std::unique_ptr<task> create_group_task(
		const std::string& input,
		const subtasks_t& subtasks,
		const_iter& a_subtask_group )
{
	std::string subgroup_type_id = *a_subtask_group;
	auto next_subtask = a_subtask_group;
	SPLOG("composite_task::create_group_task> for: %s", subgroup_type_id.c_str());
	std::vector<std::string> all_subtasks;

	while( (++next_subtask != subtasks.cend()) && (*next_subtask != subgroup_type_id) )
	{
		SPLOG_P(LOG_DEBUG, "composite_task::create_group_task> adding: %s", (*next_subtask).c_str());
		all_subtasks.push_back(*next_subtask);
	}


	// the following 2 throws are assertions.
	// If they happen, that means something is logically inconsistent in this class
	// regardless of its input.
	if ( next_subtask == subtasks.cend() )
	{
		throw io_exception("invalid subtask group");
	}
	if ( *next_subtask != subgroup_type_id )
	{
		throw io_exception("invalid subtask group: expected " + subgroup_type_id + " got: " + (*next_subtask));
	}

	// this allows the parent task to jump to the next subtask once this one is created
	a_subtask_group = next_subtask;

	if( starts_with("parallel", *a_subtask_group) )
	{
		std::unique_ptr<composite_task<true>> new_subtask = make_unique<composite_task<true>>(all_subtasks);
		new_subtask->input = input;
		return std::move(new_subtask);
	}
	else if ( starts_with("serial", *a_subtask_group) )
	{
		std::unique_ptr<composite_task<false>> new_subtask = make_unique<composite_task<false>>(all_subtasks);
		new_subtask->input = input;
		return std::move(new_subtask);
	}
	else
	{
		throw io_exception("invalid task type");
	}
}

auto do_nothing = []( const_iter& ) {};

size_t count_children( const subtasks_t& subtasks )
{
	size_t count(0LU);

	auto increase_count = [&]( const_iter& ) { count++; };

	for_each_child_in( subtasks, increase_count, do_nothing, increase_count );

	return count;
};
