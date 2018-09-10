"""Sample 2: loop compare"""
import json
import matplotlib.pyplot as plt
import poropyck

# read input from JSON
with open('sample2_input.json') as jsonfile:
    data = json.load(jsonfile)

# put files into (template, query) tuples for DTW
templates = data['waves'][:-1]  # last file is never used as template
queries = data['waves'][1:]  # first file is never used as query
picks = zip(templates, queries)

# loop over each (template, query) tuple
data['picks'] = []
for template, query in picks:
    # run poropyck
    dtw = poropyck.DTW(template['file'], query['file'], data['lengths'])
    data['picks'].append(dtw.pick())

# write output to JSON
with open('sample2_output.json', 'w') as jsonfile:
    json.dump(data, jsonfile, indent=2)

# plot velocities with error
plt.errorbar(
    [pick['velocity'] for pick in data['picks']],
    [w['pressure'] for w in data['waves'][1:]],
    xerr=[pick['velocity_error'] for pick in data['picks']]
)
plt.show()
