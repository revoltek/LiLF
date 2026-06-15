--[[
 This is a LBA AOFlagger strategy for combined subbands, version 2021-03-30
 Author: André Offringa
]]

aoflagger.require_min_version("3.1")

function execute(input)
  --
  -- Generic settings
  --

  -- What polarizations to flag? Default: input:get_polarizations() (=all that are in the input data)
  -- Other options are e.g.:
  -- { 'XY', 'YX' } to flag only XY and YX, or
  -- { 'I', 'Q' } to flag only on Stokes I and Q
  local flag_polarizations = input:get_polarizations()

  local base_threshold = 1.2 -- lower means more sensitive detection
  -- How to flag complex values, options are: phase, amplitude, real, imaginary, complex
  -- May have multiple values to perform detection multiple times
  local flag_representations = { "amplitude" }
  local iteration_count = 5 -- how many iterations to perform?
  local threshold_factor_step = 2.0 -- How much to increase the sensitivity each iteration?
  -- If the following variable is true, the strategy will consider existing flags
  -- as bad data. It will exclude flagged data from detection, and make sure that any existing
  -- flags on input will be flagged on output. If set to false, existing flags are ignored.
  local use_input_flags = true
  local transient_threshold_factor = 1.0 -- decreasing this value makes detection of transient RFI more aggressive

  --
  -- End of generic settings
  --

  local inpPolarizations = input:get_polarizations()

  if not use_input_flags then
    input:clear_mask()
  end
  -- For collecting statistics. Note that this is done after clear_mask(),
  -- so that the statistics ignore any flags in the input data.
  local copy_of_input = input:copy()

  aoflagger.normalize_bandpass(input)

  for ipol, polarization in ipairs(flag_polarizations) do
    local pol_data = input:convert_to_polarization(polarization)
    local converted_data
    local converted_copy

    for _, representation in ipairs(flag_representations) do
      converted_data = pol_data:convert_to_complex(representation)
      converted_copy = converted_data:copy()

      for i = 1, iteration_count - 1 do
        local threshold_factor = threshold_factor_step ^ (iteration_count - i)

        local sumthr_level = threshold_factor * base_threshold
        if use_input_flags then
          aoflagger.sumthreshold_masked(
            converted_data,
            converted_copy,
            sumthr_level,
            sumthr_level * transient_threshold_factor,
            true,
            true
          )
        else
          aoflagger.sumthreshold(converted_data, sumthr_level, sumthr_level * transient_threshold_factor, true, true)
        end

        -- Do timestep & channel flagging
        local chdata = converted_data:copy()
        aoflagger.threshold_timestep_rms(converted_data, 3.5)
        aoflagger.threshold_channel_rms(chdata, 3.0 * threshold_factor, true)
        converted_data:join_mask(chdata)

        -- High pass filtering steps
        converted_data:set_visibilities(converted_copy)
        if use_input_flags then
          converted_data:join_mask(converted_copy)
        end

        aoflagger.low_pass_filter(converted_data, 21, 31, 1.0, 1.0)

        -- In case this script is run from inside rfigui, calling
        -- the following visualize function will add the current result
        -- to the list of displayable visualizations.
        -- If the script is not running inside rfigui, the call is ignored.
        aoflagger.visualize(converted_data, "Fit #" .. i, i - 1)

        local tmp = converted_copy - converted_data
        tmp:set_mask(converted_data)
        converted_data = tmp

        aoflagger.visualize(converted_data, "Residual #" .. i, i + iteration_count)
        aoflagger.set_progress((ipol - 1) * iteration_count + i, #flag_polarizations * iteration_count)
      end -- end of iterations

      if use_input_flags then
        aoflagger.sumthreshold_masked(
          converted_data,
          converted_copy,
          base_threshold,
          base_threshold * transient_threshold_factor,
          true,
          true
        )
      else
        aoflagger.sumthreshold(converted_data, base_threshold, base_threshold * transient_threshold_factor, true, true)
      end
    end -- end of complex representation iteration

    if use_input_flags then
      converted_data:join_mask(converted_copy)
    end

    -- Helper function used below
    function contains(arr, val)
      for _, v in ipairs(arr) do
        if v == val then
          return true
        end
      end
      return false
    end

    if contains(inpPolarizations, polarization) then
      if input:is_complex() then
        converted_data = converted_data:convert_to_complex("complex")
      end
      input:set_polarization_data(polarization, converted_data)
    else
      input:join_mask(converted_data)
    end

    aoflagger.visualize(converted_data, "Residual #" .. iteration_count, 2 * iteration_count)
    aoflagger.set_progress(ipol, #flag_polarizations)
  end -- end of polarization iterations

  if use_input_flags then
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
