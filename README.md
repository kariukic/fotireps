## FOTIREPS - an RTS EoR Processing Suite, for Friends Of The Ionosphere. 
        __meant to save us from near 40 repetitions of trying to run all the steps__

Fotireps is a ***user-friendly*** conglomerate of different data processing tools used by the `MWA EoR group`. It was originally started as a tool for ionospheric analysis, hence the name.

Fotireps is hoped to be ***highly modularized and robust*** with each step being runnable in isolation but can also perform ***standard end to end EoR analysis*** on MWA observations.

The key functionalities are described below. A \* indicates tasks that are either still under development or not yet integrated.

### 1. Data selection* 
- Custom made functionality based on `MWA ASVO` and `MWA web services`.

### 2. Data downloading*

- Using `Giant Squid`

### 3. Flagging and flags statistics.*

- Currently done suing `Cotter`/`AOFlagger`.

### 4. DI (patch) and DD (peel) calibration.
- This is done using the `RTS`, `sourcelist_by_beam`, ...

### 5. Ionospheric analysis 

- This was the original purpose for Fotireps using tools such as `CTHULHU`,  `Offset_offsets/OI`, ...
- Requires ionospheric information from the RTS calibration

### 5. Power spectrum (PS) estimation

- Uses `CHIPS`.
- Includes power spectrum plots for individual or multiple observations.


## Dependencies:
