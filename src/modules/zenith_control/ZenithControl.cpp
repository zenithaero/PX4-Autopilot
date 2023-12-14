#include "ZenithControl.hpp"
#include "PX4Controller.h"

using namespace time_literals;
using namespace matrix;

using math::constrain;
using math::interpolate;
using math::radians;

ZenithControl::ZenithControl() :
	ModuleParams(nullptr),
	ScheduledWorkItem(MODULE_NAME, px4::wq_configurations::zenith_controller),
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
	// if (!_vehicle_angular_velocity_sub.registerCallback()) {
	// 	PX4_ERR("callback registration failed");
	// 	return false;
	// }

	// Initialize model
	PX4Controller_initialize();
	// kick off work queue
	ScheduleNow();
	return true;
}

void ZenithControl::setTrack(Vector3f currentPos, float currentYaw)
{
	_trackEnabled = true;
	_trackPos = currentPos;
	_trackYaw = currentYaw;
	PX4_INFO("Track set to: [%.2f %.2f %.2f] yaw: %.2f",
	 _trackPos(0), _trackPos(1), _trackPos(2), _trackYaw);
}


void ZenithControl::Run()
{
	if (should_exit()) {
		// _vehicle_angular_velocity_sub.unregisterCallback();
		exit_and_cleanup();
		return;
	}
	perf_begin(_loop_perf);

	// Populate input
	manual_control_setpoint_s _manual_control_setpoint{};
	if (_manual_control_setpoint_sub.copy(&_manual_control_setpoint)) {
		PX4Controller_U.CmdBusIn.rc[0] = (_manual_control_setpoint.throttle + 1.f) * 0.5f;
		PX4Controller_U.CmdBusIn.rc[1] = _manual_control_setpoint.roll;
		PX4Controller_U.CmdBusIn.rc[2] = _manual_control_setpoint.pitch;
		PX4Controller_U.CmdBusIn.rc[3] = _manual_control_setpoint.yaw;
	}

	// Populate position readings
	vehicle_local_position_s _local_pos{};
	_local_pos_sub.copy(&_local_pos);
	Vector3f currentPos(_local_pos.x, _local_pos.y, _local_pos.z);

	PX4Controller_U.EstBus_m.vNedEst[0] = _local_pos.vx;
	PX4Controller_U.EstBus_m.vNedEst[1] = _local_pos.vy;
	PX4Controller_U.EstBus_m.vNedEst[2] = _local_pos.vz;


	PX4Controller_U.EstBus_m.xNedEst[0] = _local_pos.x;
	PX4Controller_U.EstBus_m.xNedEst[1] = _local_pos.y;
	PX4Controller_U.EstBus_m.xNedEst[2] = _local_pos.z;

	// Height & height rate

	airspeed_validated_s airspeed_validated{};
	_airspeed_validated_sub.copy(&airspeed_validated);


	float tas = airspeed_validated.true_airspeed_m_s;
	PX4Controller_U.EstBus_m.hDotEst = -_local_pos.vz;
	PX4Controller_U.EstBus_m.hEst = -_local_pos.z;
	PX4Controller_U.EstBus_m.tasEst = tas;
	PX4Controller_U.EstBus_m.tasDotEst = 0;
	PX4Controller_U.EstBus_m.alphaEst = 0;
	PX4Controller_U.EstBus_m.betaEst = 0;

	// Populate track commands
	PX4Controller_U.CmdBusIn.tasCmd = _trackTas;
	PX4Controller_U.CmdBusIn.tasDotCmd = 0;
	PX4Controller_U.CmdBusIn.hCmd = -_trackPos(2);
	PX4Controller_U.CmdBusIn.hDotCmd = 0;
	PX4Controller_U.CmdBusIn.trackYaw = _trackYaw;
	PX4Controller_U.CmdBusIn.trackNE[0] = _trackPos(0);
	PX4Controller_U.CmdBusIn.trackNE[1] = _trackPos(1);


	// Position command
	// position_setpoint_triplet_s pos_sp_triplet{};
	// _pos_sp_triplet_sub.copy(&pos_sp_triplet);
	// bool prevValid = PX4_ISFINITE(pos_sp_triplet.previous.lat)
	// 		&& PX4_ISFINITE(pos_sp_triplet.previous.lon)
	// 		&& PX4_ISFINITE(pos_sp_triplet.previous.alt);
	// bool currValid = PX4_ISFINITE(pos_sp_triplet.current.lat)
	// 		&& PX4_ISFINITE(pos_sp_triplet.current.lon)
	// 		&& PX4_ISFINITE(pos_sp_triplet.current.alt);

	// if (pos_sp_triplet.previous.valid && pos_sp_triplet.current.valid) {
	// }


	//  PX4_INFO("hCmd: %.2f hEst %.2f tasEst %.2f", hCmd, -_local_pos.z, tas);


	// float current_altitude = -_local_pos.z + _local_pos.ref_alt; // Altitude AMSL in meters



	// Populate attitude readings
	vehicle_attitude_s att{};
	_att_sub.copy(&att);
	_R = matrix::Quatf(att.q);
	const matrix::Eulerf euler_angles(_R);
	PX4Controller_U.EstBus_m.eulerBodyEst[0] = euler_angles.phi();
	PX4Controller_U.EstBus_m.eulerBodyEst[1] = euler_angles.theta();
	PX4Controller_U.EstBus_m.eulerBodyEst[2] = euler_angles.psi();


	//  PX4_INFO("x y z %.2f %.2f %.2f yaw: %.2f", _local_pos.x, _local_pos.y, _local_pos.z, euler_angles.psi());


	// Populate rate readings
	vehicle_angular_velocity_s angular_velocity{};
	_vehicle_angular_velocity_sub.copy(&angular_velocity);
	Vector3f rates(angular_velocity.xyz);
	PX4Controller_U.EstBus_m.wBodyEst[0] = rates(0);
	PX4Controller_U.EstBus_m.wBodyEst[1] = rates(1);
	PX4Controller_U.EstBus_m.wBodyEst[2] = rates(2);



	// Reset manual mode
	PX4Controller_U.CmdBusIn.manualAttitude = false;
	PX4Controller_U.CmdBusIn.manualRate = false;
	PX4Controller_U.CmdBusIn.manualFM = false;
	PX4Controller_U.CmdBusIn.manualActuation = false;

	vehicle_control_mode_s	_vcontrol_mode{};
	_vehicle_control_mode_sub.copy(&_vcontrol_mode);
	// Manual pass-through
	if (_vcontrol_mode.flag_control_position_enabled) {
		// Not manual
		if (!_trackEnabled) setTrack(currentPos, euler_angles.psi());

	} else {
		// Manual control mode
		_trackEnabled = false;
		if (_vcontrol_mode.flag_control_attitude_enabled) {
			// PX4_INFO("Manual attitude control");
			PX4Controller_U.CmdBusIn.manualAttitude = true;
		} else if (_vcontrol_mode.flag_control_rates_enabled) {
			// PX4_INFO("Manual rate control");
			PX4Controller_U.CmdBusIn.manualRate = true;
		}
		else {
			// PX4_INFO("Manual control");
			PX4Controller_U.CmdBusIn.manualActuation = true;
		}
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

	// static int counter = 0;
	// if (counter ++ % 100 == 0) {
	// 	PX4_INFO("Actuation [%.2f %.2f %.2f %.2f %.2f %.2f]",
	// 		actuator_servos.control[leftAilIdx],
	// 		actuator_servos.control[rightAilIdx],
	// 		actuator_servos.control[leftFlapIdx],
	// 		actuator_servos.control[rightFlapIdx],
	// 		actuator_servos.control[elevIdx],
	// 		actuator_servos.control[rudIdx]);
	// }



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
		PX4_ERR("alloc failed!");
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
