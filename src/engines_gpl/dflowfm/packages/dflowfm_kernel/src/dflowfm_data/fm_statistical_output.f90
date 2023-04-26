module fm_statistical_output
   use m_output_config
   use m_statistical_output
   implicit none
   
   private
   type(t_output_quantity_config), allocatable, public, dimension(:) :: out_quan_conf_his
   type(t_output_quantity_config), allocatable, public, dimension(:) :: out_quan_conf_map
   type(t_output_quantity_config), allocatable, public, dimension(:) :: out_quan_conf_classmap

   type(t_output_variable_set), allocatable, public, dimension(:) :: out_variable_set_his
   type(t_output_variable_set), allocatable, public, dimension(:) :: out_variable_set_map
   type(t_output_variable_set), allocatable, public, dimension(:) :: out_variable_set_classmap
   
   public default_fm_statistical_output

   interface add_additional_entry
      module procedure add_additional_entry_strings
      module procedure add_additional_entry_doubles
      module procedure add_additional_entry_integers
      module procedure add_additional_entry_string
      module procedure add_additional_entry_double
      module procedure add_additional_entry_integer
   end interface

contains
   
   subroutine default_fm_statistical_output()
      use coordinate_reference_system
      use netcdf_utils

      integer numhis, nummap, numclass
      allocate(out_quan_conf_his(1000))
      allocate(out_quan_conf_map(1000))
      allocate(out_quan_conf_classmap(1000))

      numhis = 0
      nummap = 0
      numclass = 0

      call setoutval(out_quan_conf_his, numhis, IDX_HIS_VOLTOT,                                             'Wrihis_balance',       'total_volume',              '', '', 'TJ', UNC_LOC_WB)                                                 
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_STOR,                                               'Wrihis_balance',       'storage',                   '', '', 'TJ', UNC_LOC_WB)                
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_VOLERR,                                             'Wrihis_balance',       'volume_error',              '', '', 'TJ', UNC_LOC_WB)                                                 
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_BNDIN,                                              'Wrihis_balance',       'boundaries_in',             '', '', 'TJ', UNC_LOC_WB)                                                 
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_BNDOUT,                                             'Wrihis_balance',       'boundaries_out',            '', '', 'TJ', UNC_LOC_WB)                                                 
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_BNDTOT,                                             'Wrihis_balance',       'boundaries_total',          '', '', 'TJ', UNC_LOC_WB)                                                    
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_EXCHIN,                                             'Wrihis_balance',       'exchange_with_1D_in',       '', '', 'TJ', UNC_LOC_WB)                                                       
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_EXCHOUT,                                            'Wrihis_balance',       'exchange_with_1D_out',      '', '', 'TJ', UNC_LOC_WB)                                                       
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_EXCHTOT,                                            'Wrihis_balance',       'exchange_with_1D_total',    '', '', 'TJ', UNC_LOC_WB)                                                          
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_PRECIP_TOTAL,                                       'Wrihis_balance',       'precipitation_total',       '', '', 'TJ', UNC_LOC_WB)                                                       
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_EVAP,                                               'Wrihis_balance',       'evaporation',               '', '', 'TJ', UNC_LOC_WB)                                  
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_SOUR,                                               'Wrihis_balance',       'source_sink',               '', '', 'TJ', UNC_LOC_WB)                                        
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_InternalTidesDissipation,                           'Wrihis_balance',       'InternalTidesDissipation',  '', '', 'TJ', UNC_LOC_WB)                                                             
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_GravInput,                                          'Wrihis_balance',       'Gravitational_Input',       '', '', 'TJ', UNC_LOC_WB)                                                       
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_SalInput,                                           'Wrihis_balance',       'SAL_Input',                 '', '', 'TJ', UNC_LOC_WB)                                  
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_SalInput2,                                          'Wrihis_balance',       'SAL_Input_2''', '', 'TJ', UNC_LOC_WB)
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_GRWIN,                                              'Wrihis_balance',       'groundwater_in',            '', '', 'TJ', UNC_LOC_WB)                                                 
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_GRWOUT,                                             'Wrihis_balance',       'groundwater_out',           '', '', 'TJ', UNC_LOC_WB)                                                    
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_GRWTOT,                                             'Wrihis_balance',       'groundwater_total',         '', '', 'TJ', UNC_LOC_WB)                                                    
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_LATIN,                                              'Wrihis_balance',       'laterals_in''', '', 'TJ', UNC_LOC_WB)
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_LATOUT,                                             'Wrihis_balance',       'laterals_out',              '', '', 'TJ', UNC_LOC_WB)                                                 
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_LATTOT,                                             'Wrihis_balance',       'laterals_total',            '', '', 'TJ', UNC_LOC_WB)                                                 
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_LATIN1D,                                            'Wrihis_balance',       'laterals_in_1D',            '', '', 'TJ', UNC_LOC_WB)                                                 
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_LATOUT1D,                                           'Wrihis_balance',       'laterals_out_1D',           '', '', 'TJ', UNC_LOC_WB)                                                    
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_LATTOT1D,                                           'Wrihis_balance',       'laterals_total_1D',         '', '', 'TJ', UNC_LOC_WB)                                                    
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_LATIN2D,                                            'Wrihis_balance',       'laterals_in_2D',            '', '', 'TJ', UNC_LOC_WB)                                                 
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_LATOUT2D,                                           'Wrihis_balance',       'laterals_out_2D',           '', '', 'TJ', UNC_LOC_WB)                                                    
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_LATTOT2D,                                           'Wrihis_balance',       'laterals_total_2D',         '', '', 'TJ', UNC_LOC_WB)                                                    
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_EXTIN,                                              'Wrihis_balance',       'Qext_in',                   '', '', 'TJ', UNC_LOC_WB)                                     
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_EXTOUT,                                             'Wrihis_balance',       'Qext_out',                  '', '', 'TJ', UNC_LOC_WB)                                        
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_EXTTOT,                                             'Wrihis_balance',       'Qext_total',                '', '', 'TJ', UNC_LOC_WB)                                              
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_EXTIN1D,                                            'Wrihis_balance',       'Qext_in_1D',                '', '', 'TJ', UNC_LOC_WB)                                        
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_EXTOUT1D,                                           'Wrihis_balance',       'Qext_out_1D',               '', '', 'TJ', UNC_LOC_WB)                                     
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_EXTTOT1D,                                           'Wrihis_balance',       'Qext_total_1D',             '', '', 'TJ', UNC_LOC_WB)                                                 
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_EXTIN2D,                                            'Wrihis_balance',       'Qext_in_2D',                '', '', 'TJ', UNC_LOC_WB)                                                    
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_EXTOUT2D,                                           'Wrihis_balance',       'Qext_out_2D',               '', '', 'TJ', UNC_LOC_WB)                                           
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_EXTTOT2D,                                           'Wrihis_balance',       'Qext_total_2D',             '', '', 'TJ', UNC_LOC_WB)                                                 
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_ICEPT,                                              'Wrihis_balance',       'total_volume_interception', '', '', 'TJ', UNC_LOC_WB)                                                             
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_EVAP_ICEPT,                                         'Wrihis_balance',       'evaporation_interception',  '', '', 'TJ', UNC_LOC_WB)                                                             
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_PRECIP_GROUND,                                      'Wrihis_balance',       'precipitation_on_ground',   '', '', 'TJ', UNC_LOC_WB)                                                          
  
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_SOURCE_SINK_PRESCRIBED_DISCHARGE,                   'Wrihis_sourcesink',    'source_sink_prescribed_discharge',                   '',                                                      '',                     'm3 s-1',   UNC_LOC_SOSI   )
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_SOURCE_SINK_PRESCRIBED_SALINITY_INCREMENT,          'Wrihis_sourcesink',    'source_sink_prescribed_salinity_increment',          '',                                                      '',                     '1e-3',     UNC_LOC_SOSI   )
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_SOURCE_SINK_PRESCRIBED_TEMPERATURE_INCREMENT,       'Wrihis_sourcesink',    'source_sink_prescribed_temperature_increment',       '',                                                      '',                     'degC',     UNC_LOC_SOSI   )
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_SOURCE_SINK_CURRENT_DISCHARGE,                      'Wrihis_sourcesink',    'source_sink_current_discharge',                      '',                                                      '',                     'm3 s-1',   UNC_LOC_SOSI   )
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_SOURCE_SINK_CUMULATIVE_VOLUME,                      'Wrihis_sourcesink',    'source_sink_cumulative_volume',                      '',                                                      '',                     'm3',       UNC_LOC_SOSI   )
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_SOURCE_SINK_DISCHARGE_AVERAGE ,                     'Wrihis_sourcesink',    'source_sink_discharge_average' ,                     '',                                                      '',                     'm3 s-1',   UNC_LOC_SOSI   )
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_GENERAL_STRUCTURE_DISCHARGE,                        'Wrihis_structure_gen', 'general_structure_discharge',                        'Total discharge through general structure',             '',                     'm3 s-1',   UNC_LOC_GENSTRU)
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_GENERAL_STRUCTURE_CREST_LEVEL,                      'Wrihis_structure_gen', 'general_structure_crest_level',                      'Crest level of general structure',                      '',                     'm',        UNC_LOC_GENSTRU)
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_GENERAL_STRUCTURE_GATE_LOWER_EDGE_LEVEL,            'Wrihis_structure_gen', 'general_structure_gate_lower_edge_level',            'Gate lower edge level of general structure',            '',                     'm',        UNC_LOC_GENSTRU)
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_GENERAL_STRUCTURE_GATE_OPENING_WIDTH,               'Wrihis_structure_gen', 'general_structure_gate_opening_width',               'Gate opening width of general structure',               '',                     'm',        UNC_LOC_GENSTRU)
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_GENERAL_STRUCTURE_S1UP,                             'Wrihis_structure_gen', 'general_structure_s1up',                             'Water level upstream of general structure',             'sea_surface_height',   'm',        UNC_LOC_GENSTRU)
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_GENERAL_STRUCTURE_S1DN,                             'Wrihis_structure_gen', 'general_structure_s1dn',                             'Water level downstream of general structure',           'sea_surface_height',   'm',        UNC_LOC_GENSTRU)
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_GENERAL_STRUCTURE_HEAD,                             'Wrihis_structure_gen', 'general_structure_head',                             'Head difference across general structure',              '',                     'm',        UNC_LOC_GENSTRU)
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_GENERAL_STRUCTURE_FLOW_AREA,                        'Wrihis_structure_gen', 'general_structure_flow_area',                        'Flow area at general structure',                        '',                     'm2',       UNC_LOC_GENSTRU)
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_GENERAL_STRUCTURE_VELOCITY,                         'Wrihis_structure_gen', 'general_structure_velocity',                         'Velocity through general structure',                    '',                     'm s-1',    UNC_LOC_GENSTRU)
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_GENERAL_STRUCTURE_CREST_WIDTH,                      'Wrihis_structure_gen', 'general_structure_crest_width',                      'Crest width of general structure',                      '',                     'm',        UNC_LOC_GENSTRU)
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_GENERAL_STRUCTURE_DISCHARGE_THROUGH_GATE_OPENING,   'Wrihis_structure_gen', 'general_structure_discharge_through_gate_opening',   'Discharge through gate opening of general structure',   '',                     'm3 s-1',   UNC_LOC_GENSTRU)
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_GENERAL_STRUCTURE_DISCHARGE_OVER_GATE,              'Wrihis_structure_gen', 'general_structure_discharge_over_gate',              'Discharge over gate of general structure',              '',                     'm3 s-1',   UNC_LOC_GENSTRU)
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_GENERAL_STRUCTURE_DISCHARGE_UNDER_GATE,             'Wrihis_structure_gen', 'general_structure_discharge_under_gate',             'Discharge under gate of general structure',             '',                     'm3 s-1',   UNC_LOC_GENSTRU)
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_GENERAL_STRUCTURE_GATE_OPENING_HEIGHT,              'Wrihis_structure_gen', 'general_structure_gate_opening_height',              'Gate opening height of general structure',              '',                     'm',        UNC_LOC_GENSTRU)
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_GENERAL_STRUCTURE_GATE_UPPER_EDGE_LEVEL,            'Wrihis_structure_gen', 'general_structure_gate_upper_edge_level',            'Gate upper edge level of general structure',            '',                     'm',        UNC_LOC_GENSTRU)
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_GENERAL_STRUCTURE_VELOCITY_THROUGH_GATE_OPENING,    'Wrihis_structure_gen', 'general_structure_velocity_through_gate_opening',    'Velocity through gate opening of general structure',    '',                     'm s-1',    UNC_LOC_GENSTRU)
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_GENERAL_STRUCTURE_VELOCITY_OVER_GATE,               'Wrihis_structure_gen', 'general_structure_velocity_over_gate',               'Velocity over gate of general structure',               '',                     'm s-1',    UNC_LOC_GENSTRU)
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_GENERAL_STRUCTURE_VELOCITY_UNDER_GATE,              'Wrihis_structure_gen', 'general_structure_velocity_under_gate',              'Flow area in gate opening of general structure',        '',                     'm s-1',    UNC_LOC_GENSTRU)
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_GENERAL_STRUCTURE_FLOW_AREA_IN_GATE_OPENING,        'Wrihis_structure_gen', 'general_structure_flow_area_in_gate_opening',        'Flow area in gate opening of general structure',        '',                     'm2',       UNC_LOC_GENSTRU)
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_GENERAL_STRUCTURE_FLOW_AREA_OVER_GATE,              'Wrihis_structure_gen', 'general_structure_flow_area_over_gate',              'Flow area over gate of general structure',              '',                     'm2',       UNC_LOC_GENSTRU)
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_GENERAL_STRUCTURE_FLOW_AREA_UNDER_GATE,             'Wrihis_structure_gen', 'general_structure_flow_area_under_gate',             'Flow area under gate of general structure',             '',                     'm2',       UNC_LOC_GENSTRU)
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_GENERAL_STRUCTURE_STATE,                            'Wrihis_structure_gen', 'general_structure_state',                            'Flow state at general structure',                       '',                     '-',        UNC_LOC_GENSTRU)
      allocate(out_quan_conf_his(numhis)%additional_attributes(4))
      out_quan_conf_his(numhis)%num_additional_attributes = 4
      call ncu_add_att(out_quan_conf_his(numhis)%additional_attributes(1), 'flag_values', (/ 0, 1, 2, 3, 4 /))
      call ncu_add_att(out_quan_conf_his(numhis)%additional_attributes(2), 'flag_meanings', 'no_flow weir_free weir_submerged gate_free gate_submerged')
      call ncu_add_att(out_quan_conf_his(numhis)%additional_attributes(3), 'valid_range', (/ 0, 4 /))
      call ncu_add_att(out_quan_conf_his(numhis)%additional_attributes(4), '_FillValue', int(dmiss))
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_GENERAL_STRUCTURE_S1_ON_CREST,                      'Wrihis_structure_gen', 'general_structure_s1_on_crest',                      'Water level on crest of general structure',             '',                     'm',        UNC_LOC_GENSTRU)
      call setoutval(out_quan_conf_his, numhis, IDX_HIS_GENERAL_STRUCTURE_FORCE_DIFFERENCE,                 'Wrihis_structure_gen', 'general_structure_force_difference',                 'Force difference per unit at general structure',        '',                     'N m-1',    UNC_LOC_GENSTRU)

      call setoutval(out_quan_conf_his, numhis, idx_his_,'Wrihis_structure_dam', jahiscdam, success)
      call setoutval(out_quan_conf_his, numhis, idx_his_, 'Wrihis_structure_dam', '',                            '',                       '',                     '-',        UNC_LOC_)

      call setoutval(out_quan_conf_his, numhis, idx_his_,'Wrihis_structure_pump', jahispump, success)
      call setoutval(out_quan_conf_his, numhis, idx_his_, 'Wrihis_structure_pump', '',                            '',                       '',                     '-',        UNC_LOC_)

      call setoutval(out_quan_conf_his, numhis, idx_his_,'Wrihis_structure_gate', jahisgate, success)
      call setoutval(out_quan_conf_his, numhis, idx_his_, '', '',                            '',                       '',                     '-',        UNC_LOC_)

      call setoutval(out_quan_conf_his, numhis, idx_his_,'Wrihis_structure_weir', jahisweir, success)
      call setoutval(out_quan_conf_his, numhis, idx_his_, '', '',                            '',                       '',                     '-',        UNC_LOC_)

      call setoutval(out_quan_conf_his, numhis, idx_his_,'Wrihis_structure_orifice', jahisorif, success)
      call setoutval(out_quan_conf_his, numhis, idx_his_, '', '',                            '',                       '',                     '-',        UNC_LOC_)

      call setoutval(out_quan_conf_his, numhis, idx_his_,'Wrihis_structure_bridge', jahisbridge, success)
      call setoutval(out_quan_conf_his, numhis, idx_his_, '', '',                            '',                       '',                     '-',        UNC_LOC_)

      call setoutval(out_quan_conf_his, numhis, idx_his_,'Wrihis_structure_culvert', jahisculv,  success)
      call setoutval(out_quan_conf_his, numhis, idx_his_, '', '',                            '',                       '',                     '-',        UNC_LOC_)

      call setoutval(out_quan_conf_his, numhis, idx_his_,'Wrihis_structure_damBreak', jahisdambreak,  success)
      call setoutval(out_quan_conf_his, numhis, idx_his_, '', '',                            '',                       '',                     '-',        UNC_LOC_)

      call setoutval( ,'Wrihis_structure_uniWeir', jahisuniweir,  success)
      call setoutval(out_quan_conf_his, numhis, idx_his_, '', '',                            '',                       '',                     '-',        UNC_LOC_)

      call setoutval( ,'Wrihis_structure_compound', jahiscmpstru,  success)
      call setoutval(out_quan_conf_his, numhis, idx_his_, '', '',                            '',                       '',                     '-',        UNC_LOC_)

      call setoutval( ,'Wrihis_structure_longculvert', jahislongculv,  success)
      call setoutval(out_quan_conf_his, numhis, idx_his_, '', '',                            '',                       '',                     '-',        UNC_LOC_)

      call setoutval( ,'Wrihis_turbulence', jahistur, success)
      call setoutval(out_quan_conf_his, numhis, idx_his_, '', '',                            '',                       '',                     '-',        UNC_LOC_)

      call setoutval( ,'Wrihis_wind', jahiswind, success)
      call setoutval(out_quan_conf_his, numhis, idx_his_, '', '',                            '',                       '',                     '-',        UNC_LOC_)

      call setoutval( ,'Wrihis_rain', jahisrain, success)
      call setoutval(out_quan_conf_his, numhis, idx_his_, '', '',                            '',                       '',                     '-',        UNC_LOC_)

      call setoutval( ,'Wrihis_infiltration', jahisinfilt, success)
      call setoutval(out_quan_conf_his, numhis, idx_his_, '', '',                            '',                       '',                     '-',        UNC_LOC_)

      call setoutval( ,'Wrihis_temperature', jahistem, success)
      call setoutval(out_quan_conf_his, numhis, idx_his_, '', '',                            '',                       '',                     '-',        UNC_LOC_)

      call setoutval( ,'Wrihis_heat_fluxes', jahisheatflux, success)
      call setoutval(out_quan_conf_his, numhis, idx_his_, '', '',                            '',                       '',                     '-',        UNC_LOC_)

      call setoutval( ,'Wrihis_heatflux', jahisheatflux, success)
      call setoutval(out_quan_conf_his, numhis, idx_his_, '', '',                            '',                       '',                     '-',        UNC_LOC_)

      call setoutval( ,'Wrihis_salinity', jahissal, success)
      call setoutval(out_quan_conf_his, numhis, idx_his_, '', '',                            '',                       '',                     '-',        UNC_LOC_)

      call setoutval(out_quan_conf_his, numhis, idx_his_, '', '',                            '',                       '',                     '-',        UNC_LOC_)
      call setoutval( ,'Wrihis_density', jahisrho, success)

      call setoutval(out_quan_conf_his, numhis, idx_his_, '', '',                            '',                       '',                     '-',        UNC_LOC_)
      call setoutval( ,'Wrihis_waterlevel_s1', jahiswatlev, success)

      call setoutval(out_quan_conf_his, numhis, idx_his_, '', '',                            '',                       '',                     '-',        UNC_LOC_)
      call setoutval( ,'Wrihis_bedlevel', jahisbedlev, success)

      call setoutval(out_quan_conf_his, numhis, idx_his_, '', '',                            '',                       '',                     '-',        UNC_LOC_)
      call setoutval( ,'Wrihis_waterdepth', jahiswatdep, success)

      call setoutval(out_quan_conf_his, numhis, idx_his_, '', '',                            '',                       '',                     '-',        UNC_LOC_)
      call setoutval( ,'Wrihis_waves', jahiswav, success)

      call setoutval( ,'Wrihis_velocity_vector', jahisvelvec, success)
      call setoutval(out_quan_conf_his, numhis, idx_his_, '', '',                            '',                       '',                     '-',        UNC_LOC_)

      call setoutval( ,'Wrihis_upward_velocity_component', jahisww, success)
      call setoutval(out_quan_conf_his, numhis, idx_his_, '', '',                            '',                       '',                     '-',        UNC_LOC_)

      call setoutval( ,'Wrihis_sediment', jahissed, success)
      call setoutval(out_quan_conf_his, numhis, idx_his_, '', '',                            '',                       '',                     '-',        UNC_LOC_)

      call setoutval( ,'Wrihis_constituents', jahisconst, success)
      call setoutval(out_quan_conf_his, numhis, idx_his_, '', '',                            '',                       '',                     '-',        UNC_LOC_)

      call setoutval( ,'Wrihis_zcor', jahiszcor, success)
      call setoutval(out_quan_conf_his, numhis, idx_his_, '', '',                            '',                       '',                     '-',        UNC_LOC_)

      call setoutval( ,'Wrihis_lateral', jahislateral, success)
      call setoutval(out_quan_conf_his, numhis, idx_his_, '', '',                            '',                       '',                     '-',        UNC_LOC_)

      call setoutval( ,'Wrihis_taucurrent', jahistaucurrent, success)
      call setoutval(out_quan_conf_his, numhis, idx_his_, '', '',                            '',                       '',                     '-',        UNC_LOC_)

      call setoutval( ,'Wrihis_velocity', jahisvelocity, success)
      call setoutval(out_quan_conf_his, numhis, idx_his_, '', '',                            '',                       '',                     '-',        UNC_LOC_)

      call setoutval( ,'Wrihis_discharge', jahisdischarge, success)
      call setoutval(out_quan_conf_his, numhis, idx_his_, '', '',                            '',                       '',                     '-',        UNC_LOC_)
      call setoutval( ,'Wrimap_waterlevel_s0', jamaps0, success)
      call setoutval( ,'Wrimap_waterlevel_s1', jamaps1, success)
      call setoutval( ,'Wrimap_evaporation', jamapevap, success)
      call setoutval( ,'Wrimap_volume1', jamapvol1, success)
      call setoutval( ,'Wrimap_waterdepth', jamaphs, success)
      call setoutval( ,'Wrimap_waterdepth_hu', jamaphu, success)
      call setoutval( ,'Wrimap_ancillary_variables', jamapanc, success)
      call setoutval( ,'Wrimap_flow_analysis', jamapFlowAnalysis, success)
      call setoutval( ,'Wrimap_flowarea_au', jamapau, success)
      call setoutval( ,'Wrimap_velocity_component_u1', jamapu1, success)
      call setoutval( ,'Wrimap_velocity_component_u0', jamapu0, success)
      call setoutval( ,'Wrimap_velocity_vector', jamapucvec, success)
      call setoutval( ,'Wrimap_velocity_magnitude', jamapucmag, success)
      call setoutval( ,'Wrimap_velocity_vectorq', jamapucqvec, success)
      call setoutval( ,'Wrimap_upward_velocity_component', jamapww1, success)
      call setoutval( ,'Wrimap_density_rho', jamaprho, success)
      call setoutval( ,'Wrimap_horizontal_viscosity_viu', jamapviu, success)
      call setoutval( ,'Wrimap_horizontal_diffusivity_diu', jamapdiu, success)
      call setoutval( ,'Wrimap_flow_flux_q1', jamapq1, success)
      call setoutval( ,'Wrimap_flow_flux_q1_main', jamapq1main, success)
      call setoutval( ,'Wrimap_fixed_weir_energy_loss', jamapfw, success)
      call setoutval( ,'Wrimap_spiral_flow', jamapspir, success)
      call setoutval( ,'Wrimap_numlimdt', jamapnumlimdt, success)
      call setoutval( ,'Wrimap_taucurrent', jamaptaucurrent, success)
      call setoutval( ,'Wrimap_z0', jamapz0, success)
      call setoutval( ,'Wrimap_salinity', jamapsal, success)
      call setoutval( ,'Wrimap_chezy', jamap_chezy_elements, success)
      call setoutval( ,'Wrimap_chezy_on_flow_links', jamap_chezy_links, success)
      call setoutval( ,'Wrimap_input_roughness', jamap_chezy_input, success)
      call setoutval( ,'Wrimap_temperature', jamaptem, success)
      call setoutval( ,'Wrimap_constituents', jamapconst, success)
      call setoutval( ,'Wrimap_sediment', jamapsed, success)
      call setoutval( ,'Wrimap_turbulence', jamaptur, success)
      call setoutval( ,'Wrimap_trachytopes', jamaptrachy, success)
      call setoutval( ,'Wrimap_calibration', jamapcali, success)
      call setoutval( ,'Wrimap_rain', jamaprain, success)
      call setoutval( ,'Wrimap_interception', jamapicept, success)
      call setoutval( ,'Wrimap_wind', jamapwind, success)
      call setoutval( ,'Wrimap_windstress', jamapwindstress, success)
      call setoutval( ,'Wrimap_heat_fluxes', jamapheatflux, success)
      call setoutval( ,'Wrimap_tidal_potential', jamaptidep, success)
      call setoutval( ,'Wrimap_sal_potential', jamapselfal, success)
      call setoutval( ,'Wrimap_internal_tides_dissipation', jamapIntTidesDiss, success)
      call setoutval( ,'Wrimap_nudging', jamapnudge, success)
      call setoutval( ,'Wrimap_waves',jamapwav, success)
      call setoutval( ,'Wrimap_DTcell',jamapdtcell, success)
      call setoutval( ,'Wrimap_time_water_on_ground', jamapTimeWetOnGround, success)
      call setoutval( ,'Wrimap_freeboard', jamapFreeboard, success)
      call setoutval( ,'Wrimap_waterdepth_on_ground', jamapDepthOnGround, success)
      call setoutval( ,'Wrimap_volume_on_ground', jamapVolOnGround, success)
      call setoutval( ,'Wrimap_total_net_inflow_1d2d', jamapTotalInflow1d2d, success)
      call setoutval( ,'Wrimap_total_net_inflow_lateral', jamapTotalInflowLat, success)
      call setoutval( ,'Wrimap_water_level_gradient', jamapS1Gradient, success)
      call setoutval( ,'Writek_CdWind', jatekcd, success)
      call setoutval( ,'Wrirst_bnd', jarstbnd, success)
      call setoutval( ,'Writepart_domain', japartdomain, success)
      call setoutval( ,'Wrimap_bnd', jamapbnd, success)
      call setoutval( ,'Wrimap_Qin', jamapqin, success)
      call setoutval( ,'Wrimap_every_dt', jaeverydt, success)
      call setoutval( ,'Richardsononoutput', jaRichardsononoutput, success)
      call setoutval( ,'Wrishp_crs', jashp_crs, success)
      call setoutval( ,'Wrishp_obs', jashp_obs, success)
      call setoutval( ,'Wrishp_weir', jashp_weir, success)
      call setoutval( ,'Wrishp_thd', jashp_thd, success)
      call setoutval( ,'Wrishp_gate', jashp_gate, success)
      call setoutval( ,'Wrishp_emb', jashp_emb, success)
      call setoutval( ,'Wrishp_fxw', jashp_fxw, success)
      call setoutval( ,'Wrishp_src', jashp_src, success)
      call setoutval( ,'Wrishp_pump', jashp_pump, success)
      call setoutval( ,'Wrishp_dryarea', jashp_dry, success)
      call setoutval( ,'wrishp_genstruc', jashp_genstruc, success)
      call setoutval( ,'WriteDFMinterpretedvalues', jawriteDFMinterpretedvalues, success)
      call setoutval( ,'WriteDetailedTimers', jawriteDetailedTimers, success)
      call setoutval( ,'MbaWriteCsv', jambawritecsv, success)
      call setoutval( ,'MbaLumpFromToMba', jambalumpmba, success)
      call setoutval( ,'MbaLumpBoundaries', jambalumpbnd, success)
      call setoutval( ,'MbaLumpSourceSinks', jambalumpsrc, success)
      call setoutval( ,'MbaLumpProcesses', jambalumpproc, success)
      call setoutval( ,'waveSwartDelwaq' , jawaveSwartDelwaq, success)
      call setoutval( ,'FullGridOutput', jafullgridoutput, success)
      call setoutval( ,'EulerVelocities', jaeulervel)
      call setoutval( ,'AvgWaveQuantities'    , jaavgwavquant)
      call setoutval( ,'MomentumBalance'    , jamombal, success)
      call setoutval( ,'AvgWaveWriteAll'    , jaavgwriteall, success)
      call setoutval( ,'AvgWaveWriteH'    , jaavgwriteH    , success)   ! height
      call setoutval( ,'AvgWaveWriteE'    , jaavgwriteE    , success)   ! energy
      call setoutval( ,'AvgWaveWriteR'    , jaavgwriteR    , success)   ! roller
      call setoutval( ,'AvgWaveWriteD'    , jaavgwriteD    , success)   ! dissipation breaking+roller
      call setoutval( ,'AvgWaveWriteCel'  , jaavgwriteCel  , success)   ! celerity+group velocity
      call setoutval( ,'AvgWaveWriteDir'  , jaavgwriteDir  , success)   ! wave direction
      call setoutval( ,'AvgWaveWriteU'    , jaavgwriteU    , success)   ! velocity + stokes
      call setoutval( ,'AvgWaveWriteF'    , jaavgwriteF    , success)   ! wave forces
      call setoutval( ,'AvgWaveWriteUrms' , jaavgwriteUrms , success)   ! orbital velocities
      call setoutval( ,'AvgWaveWriteS'    , jaavgwriteS    , success)   ! water level
      call setoutval( ,'AvgWaveWriteSigm' , jaavgwriteSigm , success)   ! frequency
      call setoutval( ,'Wrimap_wet_waterdepth_threshold', epswetout, success)
      call readClasses('WaterlevelClasses', map_classes_s1)
      call readClasses('WaterdepthClasses', map_classes_hs)
      call readClasses('VelocityMagnitudeClasses', map_classes_ucmag)
      call prop_get_double(md_ptr, 'output', 'VelocityDirectionClassesInterval', map_classes_ucdirstep, success)
      call prop_get_string(md_ptr, 'output', 'WaqHorAggr', md_waqhoraggr, success)
      call prop_get_string(md_ptr, 'output', 'WaqVertAggr', md_waqvertaggr, success)
      call prop_get_string (md_ptr, 'output', 'MapOutputTimeVector', md_mptfile, success)
   end subroutine default_fm_statistical_output

   subroutine setoutval(statout, numentries, key, name, long_name, standard_name, unit, location_specifier)
      type(t_output_quantity_config),  intent(inout) :: statout(:)
      integer,                         intent(inout) :: numentries                        
      character(len=Idlen),            intent(in   ) :: key
      character(len=Idlen),            intent(in   ) :: name
      character(len=Idlen),            intent(in   ) :: long_name
      character(len=Idlen),            intent(in   ) :: standard_name
      character(len=Idlen),            intent(in   ) :: unit
      integer,                         intent(in   ) :: location_specifier

      if (present(condition)) then
         if (.not. condition) then
            ! Nothing to do
            return
         endif
      endif

      numentries = numentries+1
      statout(numentries)%key                = key             
      statout(numentries)%name               = name            
      statout(numentries)%long_name          = long_name       
      statout(numentries)%standard_name      = standard_name   
      statout(numentries)%unit               = unit            
      statout(numentries)%location_specifier = location_specifier
      statout(numentries)%num_additional_attributes = 0
   end subroutine setoutval

   subroutine setnumber_additional_entries(statout, count)
      type(t_output_quantity_config), intent(inout) :: statout
      integer                         intent(in   ) :: count         
      allocate(statout%additional_attributes(count))
      statout%num_additional_attributes = count
   end subroutine setnumber_additional_entries

   call add_additional_entry(out_quan_conf_his(numhis), 'flag_values', ints=(/ 0, 1, 2, 3, 4 /))
   call add_additional_entry(out_quan_conf_his(numhis), 'flag_meanings', chars='no_flow weir_free weir_submerged gate_free gate_submerged')
   call add_additional_entry(out_quan_conf_his(numhis), 'valid_range', (/ 0, 4 /))
   call add_additional_entry(out_quan_conf_his(numhis), '_FillValue', ints=(/int(dmiss)/))

end module fm_statistical_output