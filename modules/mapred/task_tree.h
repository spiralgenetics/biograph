#pragma once

#include "modules/mapred/task.h"
#include "modules/mapred/task_composite.h"
#include "modules/io/log.h"
#include <initializer_list>

//
// task_tree is a set of 3 classes to help build trees of tasks: 'leaf_task', 'parallel' and 'serial'
//
// so you can write:
//
// typedef leaf_task<my_func0> my_leaf_task0; // etc. with my_func1, my_func2
// REGISTER_COMPONENT(my_leaf_task0)
//
// parallel* p = new parallel {
// 				my_leaf_task0(),
//				my_leaf_task1(),
//				serial { my_other_leaf_task2(), my_other_leaf_task2() }
// };
// p->input = "yo mama so big";
// add_subtask(p);
//
//
// Its design follows the well-known composite design pattern:
// the base class is a 'component'
// A 'leaf' class is a 'component' that contains nothing but a function intended to be run as a 'task'.
// A 'composite' class is a 'component' that collects other 'components', it is also a runnable 'task'.
// There are 2 kinds of 'components': parallel and serial.
// Parallel executes its children (sounds pretty bad), concurrently.
// Serial's children are run sequentially, in case it's not evident.
//
// Input:
// A component gets its input from its parent.
// If it has no parent, you must set the input field yourself.
//
// Output:
// A parallel component's output is a msgpack-serialized vector of outputs made of the outputs of its children.
// A serial component's output is the output of its last child.
// A serial component propagates the output of the current child as the input to the next child.
//
// To use it:
//
// #1: start with your own function that follows this signature:  std::string my_func( const std::string& input )
//
// 'my_func' takes an 'input' string and returns an 'output' string
// Think of both 'input' and 'output' as either msgpack or json_serialized object, although this detail
// is completely opaque to the task_tree design.
//
//
// #2: make a 'leaf':
//
// typedef leaf_task<my_func> my_leaf;
//
// #3: start composing tasks:
//
// see first example, or the unit tests: ../test/task_tree.cpp
//

struct component_task
{
	component_task(const std::string _type):type(_type) {}
	std::string type;
	std::vector<std::string> subtasks;
	std::string input;
};

template< void (*fn)(const char* input) > class leaf_task : public task_impl<leaf_task<fn>>, public component_task
{
	public:
		leaf_task():component_task(s_type()) {}

		static std::string s_type();

		std::string type() const override { return s_type(); }

		void run()
		{
			fn(input.c_str());
			this->set_output( 0L );  // TODO: implement proper input/output management, as per documented design
		}

		TRANSFER_OBJECT
		{
			VERSION(0);
			FIELD(input, TF_STRICT);
		}
};

#define REGISTER_COMPONENT(cls) template<> std::string cls::s_type() { return #cls; }; \
	static task_register_helper<cls> _register_task ## cls;

template< bool use_parallel > class composite_task : public task_impl<composite_task<use_parallel>>, public component_task
{
	public:
		composite_task():component_task(s_type()), state(0), num_children(0) {}

		composite_task( const std::vector<std::string>& _subtasks )
			:component_task(s_type()), state(0)
		{
			subtasks.insert(subtasks.end(), _subtasks.cbegin(), _subtasks.cend());
			num_children = count_children(subtasks);
		}

		composite_task( const std::initializer_list<component_task>& list )
			: component_task(s_type()), state(0), num_children(0)
		{
			for( const auto& component : list )
			{
				if( component.type == "parallel" || component.type == "serial" )
				{
					std::string group_step_id = component.type;
					group_step_id += "_";
				       	group_step_id += std::to_string((unsigned long)&component);
					group_step_id += "_";
					group_step_id += std::to_string((unsigned long)time(0));
					subtasks.push_back( group_step_id );
					subtasks.insert(subtasks.end(), component.subtasks.cbegin(), component.subtasks.cend());
					subtasks.push_back( group_step_id );
				}
				else
				{
					subtasks.push_back( component.type );
				}
				num_children++;
			}
		}

		static std::string s_type();
		std::string type() const override { return s_type(); }

		TRANSFER_OBJECT
		{
			VERSION(0);
			FIELD(subtasks, TF_STRICT);
			FIELD(state, TF_STRICT);
			FIELD(input, TF_STRICT);
			FIELD(num_children, TF_STRICT);
		}

		void run()
		{
			if( use_parallel )
			{
				in_parallel();
			}
			else
			{
				in_series();
			}
		}

		size_t state;
		size_t num_children;

	private:

		void in_parallel()
		{
			if( state == 0 )
			{
				SPLOG_P(LOG_DEBUG, "parallel::in_parallel> creating parallel tasks");
				for ( auto it = subtasks.cbegin(); it != subtasks.cend(); )
				{
					const auto& a_task = *it;
					std::unique_ptr<task> new_task;
					if( starts_with("parallel", a_task) || starts_with("serial", a_task) )
					{
						new_task = create_group_task(input, subtasks, it);
					}
					else
					{
						std::unique_ptr<task> new_subtask = task::create_task(a_task);
						component_task* tmp = dynamic_cast<component_task*>(new_subtask.get());
						// don't forget to propagate the input
						tmp->input = input;
						new_task = std::move(new_subtask);
					}
					this->add_subtask( std::move(new_task) );
					it++;
				}
				state = 1;
				SPLOG_P(LOG_DEBUG, "parallel::in_parallel> finished creating parallel tasks");
			}
			else
			{
				SPLOG_P(LOG_DEBUG, "parallel::in_parallel> all tasks in this parallel set have been run");
				this->update_progress(1.0);
				this->set_output( 0L );
			}
		};

		void in_series()
		{
			SPLOG_P(LOG_DEBUG, "serial::in_series> begin");
			if( state == subtasks.size() )
			{
				SPLOG_P(LOG_DEBUG, "serial::in_series> done with this series of tasks");
				this->update_progress(1.0);
				this->set_output( 0L );
			}
			else
			{
				SPLOG_P(LOG_DEBUG, "serial::in_series> num_children: %lu", num_children);
				double future_progress = (num_children - state - 1)/(double)num_children;
				SPLOG_P(LOG_DEBUG, "serial::in_series> future progress: %f", future_progress);
				this->split_progress(0.0, future_progress);
				std::string new_task = subtasks[state];
				SPLOG_P(LOG_DEBUG, "serial::in_series> state=%lu, add new task: %s", state, new_task.c_str());

				std::unique_ptr<task> new_subtask;
				if( starts_with("serial", new_task) || starts_with("parallel", new_task) )
				{
					auto it = subtasks.cbegin() + state;
					new_subtask = create_group_task(input, subtasks, it);
					state = it - subtasks.cbegin();
				}
				else
				{
					std::unique_ptr<task> tmp_subtask = task::create_task( new_task );
					component_task* tmp = dynamic_cast<component_task*>(tmp_subtask.get());
					tmp->input = input;
					new_subtask = std::move(tmp_subtask);
				}
				this->add_subtask( std::move(new_subtask) );
				state++;
			}
			SPLOG_P(LOG_DEBUG, "serial::in_series> end");
		};
};

typedef composite_task<true>	parallel;
typedef composite_task<false>	serial;
