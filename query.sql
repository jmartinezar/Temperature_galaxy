SELECT TOP 150
s.z,
s.class,
g.specObjid,
g.oiii_5007_flux, g.oiii_5007_flux_err,
g.oiii_4363_flux, g.oiii_4363_flux_err,
g.oiii_4959_flux, g.oiii_4959_flux_err,
g.nii_6584_flux, g.nii_6584_flux_err,
g.nii_6548_flux, g.nii_6548_flux_err,
g.oii_3729_flux, g.oii_3729_flux_err,
g.oii_3726_flux, g.oii_3726_flux_err,
g.sii_6717_flux, g.sii_6717_flux_err,
g.sii_6731_flux, g.sii_6731_flux_err
FROM galSpecLine AS g JOIN SpecObj AS s ON s.specobjid = g.specObjID
where g.oiii_5007_flux > 0
AND g.oiii_4959_flux > 0
AND g.oiii_4363_flux > 0
AND g.oiii_5007_flux > 0
AND g.oiii_4363_flux > 0
AND g.oiii_4959_flux > 0
AND g.nii_6584_flux > 0
AND g.nii_6548_flux > 0
AND g.oii_3729_flux > 0
AND g.oii_3726_flux > 0
AND g.sii_6717_flux > 0
AND g.sii_6731_flux > 0
AND (s.class = 'GALAXY' OR s.class = 'QSO')
