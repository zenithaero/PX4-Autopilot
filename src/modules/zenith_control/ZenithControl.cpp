#include "ZenithControl.hpp"
#include "PX4Controller.h"

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
	{

	_actuator_motors_pub.advertise();
	_actuator_servos_pub.advertise();
	}

ZenithControl::~ZenithControl()
{
	perf_free(_loop_perf);
}

bool ZenithControl::init()
{
	// Initialize model
	PX4Controller_initialize();
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

	// Manual pass-through
	PX4Controller_U.CmdBusIn.manualActuation = true;

	// Populate input
	if (_manual_control_setpoint_sub.copy(&_manual_control_setpoint)) {
		PX4Controller_U.CmdBusIn.rc[0] = (_manual_control_setpoint.throttle + 1.f) * 0.5f;
		PX4Controller_U.CmdBusIn.rc[1] = _manual_control_setpoint.roll;
		PX4Controller_U.CmdBusIn.rc[2] = _manual_control_setpoint.pitch;
		PX4Controller_U.CmdBusIn.rc[3] = _manual_control_setpoint.yaw;
	}


	// Step the controller
	PX4Controller_step();

	// Retrieve output and populate uORB message
	actuator_motors_s actuator_motors;
	actuator_motors.timestamp = hrt_absolute_time();

	// Set actuator
	actuator_motors.control[0] = PX4Controller_Y.ActBus_e.motors[0];
	_actuator_motors_pub.publish(actuator_motors);


	actuator_servos_s actuator_servos;
	actuator_servos.timestamp = actuator_motors.timestamp;


	// Servos
	const size_t leftAilIdx = 0;
	const size_t rightAilIdx = 1;
	const size_t leftFlapIdx = 4;
	const size_t rightFlapIdx = 5;
	const size_t elevIdx = 2;
	const size_t rudIdx = 3;

	// actuator_servos.control[0] = 0.1;
	// actuator_servos.control[1] = 0.2;
	// actuator_servos.control[2] = 0.3;
	// actuator_servos.control[3] = 0.4;
	// actuator_servos.control[4] = 0.5;
	// actuator_servos.control[5] = 0.6;

	actuator_servos.control[leftAilIdx] = PX4Controller_Y.ActBus_e.ailerons[0];
	actuator_servos.control[leftFlapIdx] = PX4Controller_Y.ActBus_e.ailerons[1];
	actuator_servos.control[rightFlapIdx] = PX4Controller_Y.ActBus_e.ailerons[4];
	actuator_servos.control[rightAilIdx] = PX4Controller_Y.ActBus_e.ailerons[5];
	actuator_servos.control[elevIdx] = PX4Controller_Y.ActBus_e.elevators[0];
	actuator_servos.control[rudIdx] = PX4Controller_Y.ActBus_e.rudders[0];
	_actuator_servos_pub.publish(actuator_servos);



	// backup schedule
	ScheduleDelayed(10_ms);
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
