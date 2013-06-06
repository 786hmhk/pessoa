##




# Add here the directory of the cudd lib.
###########################
CUDD_DIR = ../../../cudd-2.5.0
###########################

RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: shortestPath

# Tool invocations
shortestPath: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: Cross G++ Linker'
	g++ -L$(CUDD_DIR)/mtr -L$(CUDD_DIR)/epd -L$(CUDD_DIR)/util -L$(CUDD_DIR)/st -L$(CUDD_DIR)/cudd -o "shortestPath" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(C_DEPS)$(CXX_DEPS)$(CPP_DEPS)$(C_UPPER_DEPS)$(OBJS)$(C++_DEPS)$(CC_DEPS)$(EXECUTABLES) shortestPath *.dot
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:


