#pragma once

#include <drivers/drv_hrt.h>
#include <lib/mathlib/mathlib.h>
#include <lib/parameters/param.h>
#include <lib/perf/perf_counter.h>

#include <px4_platform_common/px4_config.h>
#include <px4_platform_common/defines.h>
#include <px4_platform_common/module.h>
#include <px4_platform_common/module_params.h>
#include <px4_platform_common/posix.h>
#include <px4_platform_common/tasks.h>
#include <px4_platform_common/px4_work_queue/ScheduledWorkItem.hpp>

#include <uORB/Publication.hpp>
#include <uORB/PublicationMulti.hpp>
#include <uORB/Subscription.hpp>
#include <uORB/SubscriptionMultiArray.hpp>
#include <uORB/SubscriptionCallback.hpp>
#include <uORB/topics/actuator_controls_status.h>
#include <uORB/topics/airspeed_validated.h>
#include <uORB/topics/battery_status.h>
#include <uORB/topics/control_allocator_status.h>
#include <uORB/topics/manual_control_setpoint.h>
#include <uORB/topics/normalized_unsigned_setpoint.h>
#include <uORB/topics/parameter_update.h>
#include <uORB/topics/rate_ctrl_status.h>

#include <uORB/topics/vehicle_attitude.h>
#include <uORB/topics/vehicle_angular_velocity.h>
#include <uORB/topics/vehicle_control_mode.h>
#include <uORB/topics/vehicle_land_detected.h>
#include <uORB/topics/vehicle_rates_setpoint.h>
#include <uORB/topics/vehicle_status.h>
#include <uORB/topics/vehicle_thrust_setpoint.h>
#include <uORB/topics/vehicle_torque_setpoint.h>
#include <uORB/topics/vehicle_local_position.h>
#include <uORB/topics/vehicle_local_position_setpoint.h>
#include <uORB/topics/position_setpoint_triplet.h>
#include <uORB/topics/trajectory_setpoint.h>
#include <uORB/topics/airspeed_validated.h>

#include <uORB/topics/actuator_motors.h>
#include <uORB/topics/actuator_servos.h>


using uORB::SubscriptionData;

using namespace time_literals;

class ZenithControl final : public ModuleBase<ZenithControl>, public ModuleParams,
	public px4::ScheduledWorkItem
{
public:
	ZenithControl();
	~ZenithControl() override;

	/** @see ModuleBase */
	static int task_spawn(int argc, char *argv[]);

	/** @see ModuleBase */
	static int custom_command(int argc, char *argv[]);

	/** @see ModuleBase */
	static int print_usage(const char *reason = nullptr);

	bool init();
	void setTrack(matrix::Vector3f currentPos, float currentYaw);

private:
	void Run() override;

	// uORB::SubscriptionCallbackWorkItem _vehicle_angular_velocity_sub{this, ORB_ID(vehicle_angular_velocity)};
	// uORB::SubscriptionCallbackWorkItem _att_sub{this, ORB_ID(vehicle_attitude)};
	uORB::Subscription _local_pos_sub{ORB_ID(vehicle_local_position)};
	uORB::Subscription _vehicle_angular_velocity_sub{ORB_ID(vehicle_angular_velocity)};
	uORB::Subscription _att_sub{ORB_ID(vehicle_attitude)};

	// uORB::SubscriptionInterval _parameter_update_sub{ORB_ID(parameter_update), 1_s};
//
	// uORB::Subscription _battery_status_sub{ORB_ID(battery_status)};
	uORB::Subscription _manual_control_setpoint_sub{ORB_ID(manual_control_setpoint)};
	uORB::Subscription _trajectory_setpoint_sub{ORB_ID(trajectory_setpoint)};
	uORB::Subscription _airspeed_validated_sub{ORB_ID(airspeed_validated)};

	// uORB::Subscription _rates_sp_sub{ORB_ID(vehicle_rates_setpoint)};
	uORB::Subscription _vehicle_control_mode_sub{ORB_ID(vehicle_control_mode)};
	// uORB::Subscription _vehicle_land_detected_sub{ORB_ID(vehicle_land_detected)};
	// uORB::Subscription _vehicle_status_sub{ORB_ID(vehicle_status)};
	// uORB::Subscription _vehicle_rates_sub{ORB_ID(vehicle_angular_velocity)};
//
	uORB::SubscriptionMultiArray<control_allocator_status_s, 2> _control_allocator_status_subs{ORB_ID::control_allocator_status};
	uORB::Publication<actuator_controls_status_s>	_actuator_controls_status_pub;
	uORB::Publication<vehicle_torque_setpoint_s>	_vehicle_torque_setpoint_pub;
	uORB::Publication<vehicle_thrust_setpoint_s>	_vehicle_thrust_setpoint_pub;
//

	//
	uORB::Publication<actuator_motors_s>	_actuator_motors_pub{ORB_ID(actuator_motors)};
	uORB::Publication<actuator_servos_s>	_actuator_servos_pub{ORB_ID(actuator_servos)};


	// uORB::SubscriptionData<airspeed_validated_s> _airspeed_validated_sub{ORB_ID(airspeed_validated)};
//
	// uORB::Publication<vehicle_rates_setpoint_s>	_rate_sp_pub{ORB_ID(vehicle_rates_setpoint)};
	// uORB::PublicationMulti<rate_ctrl_status_s>	_rate_ctrl_status_pub{ORB_ID(rate_ctrl_status)};
	// uORB::Publication<normalized_unsigned_setpoint_s> _flaps_setpoint_pub{ORB_ID(flaps_setpoint)};
	// uORB::Publication<normalized_unsigned_setpoint_s> _spoilers_setpoint_pub{ORB_ID(spoilers_setpoint)};
//
	// vehicle_thrust_setpoint_s		_vehicle_thrust_setpoint{};
	// vehicle_torque_setpoint_s		_vehicle_torque_setpoint{};
	// vehicle_rates_setpoint_s		_rates_sp{};
	// vehicle_status_s			_vehicle_status{};

	perf_counter_t _loop_perf;
	matrix::Dcmf _R{matrix::eye<float, 3>()};

	// Waypoint mode
	bool _trackEnabled = false;
	matrix::Vector3f _trackPos;
	float _trackYaw = 0;
	float _trackTas = 15;
	// hrt_abstime _last_run{0};
};
