################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../ShortestPath.cpp \
../main.cpp 

CC_SRCS += \
../cuddObj.cc 

CPP_DEPS += \
./ShortestPath.d \
./main.d 

OBJS += \
./ShortestPath.o \
./cuddObj.o \
./main.o 

CC_DEPS += \
./cuddObj.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I$(CUDD_DIR)/include -O0 -g3 -Wall -c -fmessage-length=0 -mtune=native -malign-double -DHAVE_IEEE_754 -DBSD -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

%.o: ../%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I$(CUDD_DIR)/include -O0 -g3 -Wall -c -fmessage-length=0 -mtune=native -malign-double -DHAVE_IEEE_754 -DBSD -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


