px4_add_module(
	MODULE modules__zenith_control
	MAIN zenith_control
	STACK_MAX 8192
	STACK_MAIN 8192
	INCLUDES
		PX4Controller/px4/PX4Controller_grt_rtw
		PX4Controller/MATLAB_R2022a.app/extern/include
		PX4Controller/MATLAB_R2022a.app/simulink/include
		PX4Controller/MATLAB_R2022a.app/rtw/c/src
		PX4Controller/MATLAB_R2022a.app/toolbox/shared/spc/src_ml/extern/include
		PX4Controller/MATLAB_R2022a.app/bin/maci64
	SRCS
		PX4Controller/px4/PX4Controller_grt_rtw/PX4Controller.cpp
		PX4Controller/px4/PX4Controller_grt_rtw/PX4Controller_data.cpp
		PX4Controller/px4/PX4Controller_grt_rtw/rt_nonfinite.cpp
		PX4Controller/px4/PX4Controller_grt_rtw/rtGetInf.cpp
		PX4Controller/px4/PX4Controller_grt_rtw/rtGetNaN.cpp
		ZenithControl.cpp
		ZenithControl.hpp
	DEPENDS
		px4_work_queue
	COMPILE_FLAGS
		-Wno-float-equal -Wno-double-promotion
	)
