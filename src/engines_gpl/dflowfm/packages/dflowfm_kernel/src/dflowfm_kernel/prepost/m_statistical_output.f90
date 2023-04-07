module m_statistical_output
   
   contains
   
subroutine update_moving_average(stat_item)

type(t_statistical_output_item), intent(inout) :: stat_item
integer :: inew, iold

inew = stat_item%current_step
iold = MOD(stat_item%current_step+1,stat_item%total_steps_count)

if (total_steps_count >= windowsize) then !only remove oldest sample if window is full
   stat_item%moving_average = stat_item%moving_average - samples(:,iold)*timesteps(iold) + samples(:,inew)*timesteps(inew)
   stat_item%timestepsum = stat_item%stimestep_sum - timesteps(iold) + timesteps(inew)
else
   stat_item%moving_average = stat_item%moving_average + samples(:,inew)*timesteps(inew)
   stat_item%timestepsum = stat_item%stimestep_sum + timesteps(inew)
endif

end subroutine update_moving_average

subroutine update_statistical_output(stat_item)

if (stat_item%operation_id >= 2) then ! max/min of moving average requested
   call update_moving_average(stat_item)
endif

select case (stat_item%operation_id)
case (1) !SO_CURRENT
   stat_item%output = stat_item%input
case (2) !SO_AVERAGE
   stat_item%output = stat_item%output + stat_item%input * stat_item%timesteps(stat_item%current_step)
   stat_item%timestepsum = stat_item%timestepsum + stat_item%timesteps(stat_item%current_step)
case (3) !SO_MAX
   stat_item%output = max(stat_item%output,stat_item%moving_average/stat_item%timestepsum)
case (4) !SO_MIN
   stat_item%output = min(stat_item%output,stat_item%moving_average/stat_item%timestepsum)
end select

end subroutine update_statistical_output

subroutine write_statistical_output(stat_item)

if (stat_item%operation_id == 2)
   stat_item%output = stat_item%output/stat_item%timestepsum
endif

!call netcdf_write_his/map/else (stat_output_item%output)

end subroutine update_statistical_output
   
end module
   