import json

def abort(parser, msg=""):
	parser.error(msg)

def read_json(path):
	with open(path) as fid:
		return json.load(fid)

def update_list(li, indexes):
	up = 0
	for ix in indexes:
		del li[ix-up]
		up += 1
	return li
