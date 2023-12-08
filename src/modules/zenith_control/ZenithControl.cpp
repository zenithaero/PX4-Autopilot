#include "ZenithControl.hpp"

using namespace time_literals;
using namespace matrix;

using math::constrain;
using math::interpolate;
using math::radians;

ZenithControl::ZenithControl() :
	ModuleParams(nullptr),
	ScheduledWorkItem(MODULE_NAME, px4::wq_configurations::nav_and_controllers),
	_actuator_controls_status_pub(ORB_ID(actuator_controls_status_0)),
	_vehicle_torque_setpoint_pub(ORB_ID(vehicle_torque_setpoint)),
	_vehicle_thrust_setpoint_pub(ORB_ID(vehicle_thrust_setpoint)),
	_loop_perf(perf_alloc(PC_ELAPSED, MODULE_NAME": cycle"))
	{}

ZenithControl::~ZenithControl()
{
	perf_free(_loop_perf);
}

bool ZenithControl::init()
{
	// kick off work queue
	ScheduleNow();
	return true;
}


void ZenithControl::Run()
{
	if (should_exit()) {
		// _vehicle_angular_velocity_sub.unregisterCallback();
		exit_and_cleanup();
		return;
	}
	perf_begin(_loop_perf);

	// backup schedule
	ScheduleDelayed(20_ms);
	perf_end(_loop_perf);
}


int ZenithControl::task_spawn(int argc, char *argv[])
{
	ZenithControl *instance = new ZenithControl();
	if (instance) {
		_object.store(instance);
		_task_id = task_id_is_work_queue;

		if (instance->init()) {
			return PX4_OK;
		}
	} else {
		PX4_ERR("alloc failed");
	}

	delete instance;
	_object.store(nullptr);
	_task_id = -1;

	return PX4_ERROR;
}

int ZenithControl::custom_command(int argc, char *argv[])
{
	return print_usage("unknown command");
}

int ZenithControl::print_usage(const char *reason)
{
	if (reason) {
		PX4_WARN("%s\n", reason);
	}

	PRINT_MODULE_DESCRIPTION(
		R"DESCR_STR(
### Description
zenith_control is a fixed-wing controller.

)DESCR_STR");

	PRINT_MODULE_USAGE_NAME("zenith_control", "controller");
	PRINT_MODULE_USAGE_COMMAND("start");
	PRINT_MODULE_USAGE_DEFAULT_COMMANDS();

	return 0;
}

extern "C" __EXPORT int zenith_control_main(int argc, char *argv[])
{
	return ZenithControl::main(argc, argv);
}
