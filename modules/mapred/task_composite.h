#pragma once

#include "modules/mapred/task.h" 
#include <string>
#include <vector> 

//
// Just a list of helper functions to manage subtasks
//

// returns true if 'source' starts with 'prefix'
bool starts_with(const std::string& prefix, const std::string& source);

typedef std::vector<std::string> subtasks_t;
typedef subtasks_t::const_iterator const_iter;
typedef std::function< void ( const_iter& child ) > component_handler_t;

void for_each_child_in(
		const subtasks_t& subtasks,
		component_handler_t on_composite,
		component_handler_t on_composite_subtask,
		component_handler_t on_leaf );

// Creates a task with 'input' for the subgroup starting at 'a_subtask_group' in 'subtasks'
// Also modifies a_subtask_group to point to the last subtask of 'a_subtask_group'
//
std::unique_ptr<task> create_group_task(
		const std::string& input,
		const subtasks_t& subtasks,
		const_iter& a_subtask_group );

// Returns the number of children in subtasks. If children have children, those are not counted.
size_t count_children( const subtasks_t& subtasks );
