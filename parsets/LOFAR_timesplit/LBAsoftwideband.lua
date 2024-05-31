--[[
 This is a LBA AOFlagger strategy for combined subbands, version 2021-01-13
 Author: André Offringa
]]--

aoflagger.require_min_version("3.0.3")

function execute(input)

  --
  -- Generic settings
  --

  -- What polarizations to flag? Default: input:get_polarizations() (=all that are in the input data)
  -- Other options are e.g.:
  -- { 'XY', 'YX' } to flag only XY and YX, or
  -- { 'I', 'Q' } to flag only on Stokes I and Q
  local flag_polarizations = input:get_polarizations()

  local base_threshold = 1.0  -- lower means more sensitive detection
  -- How to flag complex values, options are: phase, amplitude, real, imaginary, complex
  -- May have multiple values to perform detection multiple times
  local flag_representations = { "amplitude" }
  local iteration_count = 5  -- how many iterations to perform?
  local threshold_factor_step = 2.0 -- How much to increase the sensitivity each iteration?
  -- If the following variable is true, the strategy will consider existing flags
  -- as bad data. It will exclude flagged data from detection, and make sure that any existing
  -- flags on input will be flagged on output. If set to false, existing flags are ignored.
  local exclude_original_flags = false
  local frequency_resize_factor = 1.0 -- Amount of "extra" smoothing in frequency direction
  local transient_threshold_factor = 1.0 -- decreasing this value makes detection of transient RFI more aggressive
 
  --
  -- End of generic settings
  --

  local inpPolarizations = input:get_polarizations()

  if(not exclude_original_flags) then
    input:clear_mask()
  end
  -- For collecting statistics. Note that this is done after clear_mask(),
  -- so that the statistics ignore any flags in the input data. 
  local copy_of_input = input:copy()

  aoflagger.normalize_bandpass(input)
  
  for ipol,polarization in ipairs(flag_polarizations) do
 
    local pol_data = input:convert_to_polarization(polarization)

    for _,representation in ipairs(flag_representations) do

      data = pol_data:convert_to_complex(representation)
      local original_data = data:copy()

      for i=1,iteration_count-1 do
        local threshold_factor = math.pow(threshold_factor_step, iteration_count-i)

        local sumthr_level = threshold_factor * base_threshold
        if(exclude_original_flags) then
          aoflagger.sumthreshold_masked(data, original_data, sumthr_level, sumthr_level*transient_threshold_factor, true, true)
        else
          aoflagger.sumthreshold(data, sumthr_level, sumthr_level*transient_threshold_factor, true, true)
        end

        -- Do timestep & channel flagging
        local chdata = data:copy()
        aoflagger.threshold_timestep_rms(data, 3.5)
        aoflagger.threshold_channel_rms(chdata, 3.0 * threshold_factor, true)
        data:join_mask(chdata)

        -- High pass filtering steps
        data:set_visibilities(original_data)
        if(exclude_original_flags) then
          data:join_mask(original_data)
        end

        local resized_data = aoflagger.downsample(data, 3, frequency_resize_factor, true)
        aoflagger.low_pass_filter(resized_data, 21, 31, 2.5, 5.0)
        aoflagger.upsample(resized_data, data, 3, frequency_resize_factor)

        -- In case this script is run from inside rfigui, calling
        -- the following visualize function will add the current result
        -- to the list of displayable visualizations.
        -- If the script is not running inside rfigui, the call is ignored.
        aoflagger.visualize(data, "Fit #"..i, i-1)

        local tmp = original_data - data
        tmp:set_mask(data)
        data = tmp

        aoflagger.visualize(data, "Residual #"..i, i+iteration_count)
        aoflagger.set_progress((ipol-1)*iteration_count+i, #flag_polarizations*iteration_count )
      end -- end of iterations

      if(exclude_original_flags) then
        aoflagger.sumthreshold_masked(data, original_data, base_threshold, base_threshold*transient_threshold_factor, true, true)
      else
        aoflagger.sumthreshold(data, base_threshold, base_threshold*transient_threshold_factor, true, true)
      end
    end -- end of complex representation iteration

    -- Helper function used in the strategy
    function contains(arr, val)
      for _,v in ipairs(arr) do
        if v == val then return true end
      end
      return false
    end

    if contains(inpPolarizations, polarization) then
      if input:is_complex() then
        data = data:convert_to_complex("complex")
      end
      input:set_polarization_data(polarization, data)
    else
      input:join_mask(data)
    end

    aoflagger.visualize(data, "Residual #"..iteration_count, 2*iteration_count)
    aoflagger.set_progress(ipol, #flag_polarizations )
  end -- end of polarization iterations

  if(exclude_original_flags) then
    aoflagger.scale_invariant_rank_operator_masked(input, copy_of_input, 0.2, 0.2)
  else
    aoflagger.scale_invariant_rank_operator(input, 0.2, 0.2)
  end

  aoflagger.threshold_timestep_rms(input, 4.0)

  if input:is_complex() and input:has_metadata() then
    -- This command will calculate a few statistics like flag% and stddev over
    -- time, frequency and baseline and write those to the MS. These can be
    -- visualized with aoqplot.
    aoflagger.collect_statistics(input, copy_of_input)
  end
  input:flag_nans()
end

