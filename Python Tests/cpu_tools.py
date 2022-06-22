import os
import signal
import subprocess
import time
import psutil
import numpy as np
import math

import xml.etree.cElementTree as cET
import xml.etree.ElementTree as ET
import sys

def spk_circle(N):
	layout = cET.Element("layout")
	delta = 360/N
	angle = -180.0
	for i in range(N):
    		angle += delta
    		cET.SubElement(layout, "speaker", az=f"{angle:.4f}", el="0.0000", r="1.0000")

	tree = cET.ElementTree(layout)
	tree.write("autospeakers.xml",encoding='utf-8', xml_declaration=True)
	time.sleep(0.5)


def get_plugin(plugin):

	plugin_string = "<receiver name=\"out\" type=\""+plugin+"\" layout=\"autospeakers.xml\"><position>0 0 0 0</position><orientation importcsv=\"demo_04_data/stepped_rotation.csv\" /></receiver>"
	
	
	if ("hoa3d_hybrid" in plugin) or ("hoa2d_hybrid" in plugin):
		splits = plugin.split("&")
		plugin = splits[0]
		order = splits[1]
		diviser = splits[2]
		#plugin_string = "<receiver name=\"out\" type=\""+plugin+"\" layout=\"autospeakers.xml\" order=\""+order+"\" diviser=\""+diviser+"\"><position>0 0 0 0</position><orientation importcsv=\"demo_04_data/stepped_rotation.csv\" /></receiver>"
		plugin_string = "<receiver name=\"out\" type=\""+plugin+"\" order=\""+order+"\" diviser=\""+diviser+"\"><position>0 0 0 0</position><orientation importcsv=\"demo_04_data/stepped_rotation.csv\" /></receiver>"
	elif "hoa3d_enc_lut" in plugin:
		splits = plugin.split("&")
		plugin = splits[0]
		order = splits[1]
		diviser = splits[2]
		plugin_string = "<receiver name=\"out\" type=\""+plugin+"\" order=\""+order+"\" diviser=\""+diviser+"\"><position>0 0 0 0</position><orientation importcsv=\"demo_04_data/stepped_rotation.csv\" /></receiver>"
	elif ("hoa3d_enc" in plugin) or ("hoa2d_enc" in plugin):
		splits = plugin.split("&")
		plugin = splits[0]
		order = splits[1]
		plugin_string = "<receiver name=\"out\" type=\""+plugin+"\" order=\""+order+"\"><position>0 0 0 0</position><orientation importcsv=\"demo_04_data/stepped_rotation.csv\" /></receiver>"
	elif "hoa3d" in plugin:
		splits = plugin.split("&")
		plugin = splits[0]
		order = splits[1]
		plugin_string = "<receiver name=\"out\" type=\""+plugin+"\" layout=\"autospeakers.xml\"> order=\""+order+"\"><position>0 0 0 0</position><orientation importcsv=\"demo_04_data/stepped_rotation.csv\" /></receiver>"
	
	return plugin_string
	


	
def sources_point(N,plugin):
	source_string ="<source name=\"target\"><position>0 1 0 0</position><sound><plugins><pink level=\"90\" fmax=\"16000\"/></plugins></sound></source>"
	
	plugin_string = get_plugin(plugin)
	
	et = ET.parse('scenetoedit.tsc')
	root = et.getroot()
	scene = root[0]
	plugin_el = ET.fromstring(plugin_string)
	scene.append(plugin_el)
	for i in range(N):
		position = "0 1 0 0"
		name = "target" + str(i)
		
		source_el = ET.fromstring(source_string)
		source_el.find("position").text = position
		source_el.set('name',name)

		scene.append(source_el)
	et.write('sources.tsc')
	
def sources_circle(N,plugin):
	source_string ="<source name=\"target\"><position>0 1 0 0</position><sound><plugins><pink level=\"90\" fmax=\"16000\"/></plugins></sound></source>"
	
	plugin_string = get_plugin(plugin)
	
	et = ET.parse('scenetoedit.tsc')
	root = et.getroot()
	scene = root[0]
	plugin_el = ET.fromstring(plugin_string)
	scene.append(plugin_el)
	for i in range(N):
		angle_deg = i*2*math.pi/N
		z = "0"
		x = str(math.cos(angle_deg))
		y = str(math.sin(angle_deg))
		position = "0 "+x+" "+y+" 0"
		name = "target" + str(i)
		
		source_el = ET.fromstring(source_string)
		source_el.find("position").text = position
		source_el.set('name',name)

		scene.append(source_el)
	et.write('sources.tsc')
	
	
def cpu_times(cmd,run_time):

	pro = subprocess.Popen(cmd,stdout=subprocess.PIPE, preexec_fn=os.setsid)

	ps = psutil.Process(pro.pid)
	time.sleep(run_time)
	user =ps.cpu_times().user
	system = ps.cpu_times().system
	#print(ps.cpu_times().children_user)
	os.killpg(os.getpgid(pro.pid),signal.SIGTERM)
	time.sleep(0.5)
	return user,system
	
def cmd_time(cmd):
	start = time.time()
	pro = subprocess.run(cmd,stdout=subprocess.PIPE)

	return time.time()-start

def logspace(startp,endp,n,base=10,zero=True):
	arr = list(np.logspace(startp,endp,num = n,base=base).astype(int))
	if zero:
		arr.insert(0,0)
	return arr
	
