--[[
 This is the LOFAR strategy. It is almost equal to the
 generic "minimal" AOFlagger strategy, version 2020-06-14.
 Author: André Offringa
]]--

function execute(input)

  --
  -- Generic settings
  --

  local base_threshold = 1.0  -- lower means more sensitive detection
  -- How to flag complex values, options are: phase, amplitude, real, imaginary, complex
  local representation = "amplitude"
  local iteration_count = 3  -- how many iterations to perform?
  local threshold_factor_step = 2.0 -- How much to increase the sensitivity each iteration?
  local frequency_resize_factor = 3.0 -- Amount of "extra" smoothing in frequency direction
  local transient_threshold_factor = 1.0 -- decreasing this value makes detection of transient RFI more aggressive

  --
  -- End of generic settings
  --

  local inpPolarizations = input:get_polarizations()

  input:clear_mask()

  -- For collecting statistics. Note that this is done after clear_mask(),
  -- so that the statistics ignore any flags in the input data.
  local copy_of_input = input:copy()

  for ipol,polarization in ipairs(inpPolarizations) do

    local data = input:convert_to_polarization(polarization)

    data = data:convert_to_complex(representation)
    local original_data = data:copy()

    for i=1,iteration_count-1 do
      local threshold_factor = math.pow(threshold_factor_step, iteration_count-i)

      local sumthr_level = threshold_factor * base_threshold
      aoflagger.sumthreshold(data, sumthr_level, sumthr_level*transient_threshold_factor, true, true)

      -- Do timestep & channel flagging
      local chdata = data:copy()
      aoflagger.threshold_timestep_rms(data, 3.5)
      aoflagger.threshold_channel_rms(chdata, 3.0 * threshold_factor, true)
      data:join_mask(chdata)

      -- High pass filtering steps
      data:set_visibilities(original_data)
      local resized_data = aoflagger.downsample(data, 3, frequency_resize_factor, true)
      aoflagger.low_pass_filter(resized_data, 21, 31, 2.5, 5.0)
      aoflagger.upsample(resized_data, data, 3, frequency_resize_factor)

      local tmp = original_data - data
      tmp:set_mask(data)
      data = tmp

      aoflagger.set_progress((ipol-1)*iteration_count+i, #inpPolarizations*iteration_count )
    end -- end of iterations

    aoflagger.sumthreshold(data, base_threshold, base_threshold*transient_threshold_factor, true, true)

    if input:is_complex() then
      data = data:convert_to_complex("complex")
    end
    input:set_polarization_data(polarization, data)

    aoflagger.set_progress(ipol, #inpPolarizations )
  end -- end of polarization iterations

  aoflagger.scale_invariant_rank_operator(input, 0.2, 0.2)
  aoflagger.threshold_timestep_rms(input, 4.0)
  input:flag_nans()
end

