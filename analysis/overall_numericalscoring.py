import pandas as pd
import sys

snrfile = 'numericalscoring/snr_scores.csv'
contrastfile = 'numericalscoring/contrast_scores.csv'

snrweight = float(sys.argv[1])
contrastweight = float(sys.argv[2])

fulldata = dict()
snrdf = pd.read_csv(snrfile)
contrastdf = pd.read_csv(contrastfile)

for _, row in snrdf.iterrows():
    ann, sbs, mov, spec, nb, cs, hp, score = row
    params = (ann, sbs, mov, spec, nb, cs, hp)
    fulldata[params] = list()
    fulldata[params].append(score)
for _, row in contrastdf.iterrows():
    ann, sbs, mov, spec, nb, cs, hp, score = row
    params = (ann, sbs, mov, spec, nb, cs, hp)
    fulldata[params].append(score)

for params in fulldata.keys():
    snrscore, contrastscore = fulldata[params]
    weighted_mean = (snrweight * snrscore + contrastweight * contrastscore) / (snrweight + contrastweight)
    fulldata[params] = weighted_mean

annuli, subsections, movement, spectrum, numbasis, corr_smooth, highpass, scores = list(), list(), list(), list(), \
                                                                                   list(), list(), list(), list()

for params, overall_score in zip(fulldata.keys(), fulldata.values()):
    ann, sbs, mov, spec, nb, cs, hp = params
    annuli.append(ann)
    subsections.append(sbs)
    movement.append(mov)
    spectrum.append(spec)
    numbasis.append(nb)
    corr_smooth.append(cs)
    highpass.append(hp)
    scores.append(overall_score)

overalldf = pd.DataFrame({'Annuli': annuli, 'Subsections': subsections, 'Movement': movement, 'Spectrum': spectrum,
                          'Numbasis': numbasis, 'Corr_Smooth': corr_smooth, 'Highpass': highpass, 'Score': scores})
overalldf.to_csv('numericalscoring/overall_scores.csv')