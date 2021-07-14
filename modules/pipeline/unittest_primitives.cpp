#include "modules/pipeline/primitives.h"
#include <gtest/gtest.h>

class unittest_primitives_environment : public ::testing::Environment {
	void SetUp() override {
		add_primitives();
	}
};

::testing::Environment* const unittest_primitives_env __attribute__((unused)) =
	  ::testing::AddGlobalTestEnvironment(new unittest_primitives_environment);
