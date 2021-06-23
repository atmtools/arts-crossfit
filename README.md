# arts-crossfit
ARTS HITRAN crosssection absorption model

## Calculate crosssections with the model

1. Download the precalculated coefficients files from the `arts-xml-data` package and put them in the `coefficients/` directory.

   To generate the coefficients yourself, see section ["Generate fit coefficients"](#generate-fit-coefficients) below.

2. Run the example script.
   ```
   scripts/Xsec_Calculate.py
   ```

   Creates a plot of CFC-11 crossections in `plots/Xsecs/`.

## Generate model coefficients

1. Download Hitran Crosssection data and unpack in `data/HitranXsec/`
   ```
   cd data/HitranXsec/
   curl -O https://hitran.org/data/xsec/xss.tar.gz
   tar -zxf xss.tar.gz
   cd -
   ```

2. Convert crosssection data to json format
   ```
   scripts/Xsec_ConvertHitranToJson.py
   ```

3. Generate harmonized data
   ```
   scripts/Xsec_DefineBandsAndHarmonizeData.py
   ```
   Choose option 0 to use the predefined band configuration.

4. Generate model coefficients
   ```
   scripts/Xsec_CalculateFitCoefficients.py -p
   ```

   The `-p` option generates optional diagnostic plots in `plots/SPECIES/`.

   The XML data files are stored in `coefficients/`.
