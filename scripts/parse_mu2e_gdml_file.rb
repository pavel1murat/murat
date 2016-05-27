#!/usr/bin/ruby
------------------------------------------------------------------------------
# remove address part from multiple names 
# <materials> ... </materials>
# <solids> ... </solids>
# <structure> ..</structure>
# <setup> ...</setup>
#------------------------------------------------------------------------------
require 'rexml/document'
include REXML

xmlfile = File.new("/home/murat/figures/mu2e/gdml/mu2e_geometry_v5_7_6.gdml")

# xmlfile = File.new("./test1.gdml")
doc     = Document.new(xmlfile)
root    = doc.root
#------------------------------------------------------------------------------
# 1. isotope names
#------------------------------------------------------------------------------
doc.root.elements.each("//isotope") do |e|
#  p e
  x= e.attributes["name"]
#  p "name="+x
  e.attributes["name"] = x[0,x.length-14];
#  p e.attributes["name"]
end
#------------------------------------------------------------------------------
# 2. element names
#------------------------------------------------------------------------------
doc.root.elements.each("//element") do |e|
  x= e.attributes["name"]
  e.attributes["name"] = x[0,x.length-14];
end
#------------------------------------------------------------------------------
# 3. element fraction references
#------------------------------------------------------------------------------
# doc.root.elements.each("//element/fraction") do |e|
doc.root.elements.each("//*/fraction") do |e|
  x= e.attributes["ref"]
  e.attributes["ref"] = x[0,x.length-14];
end
#------------------------------------------------------------------------------
# 4. material names
#------------------------------------------------------------------------------
doc.root.elements.each("//material") do |e|
  x= e.attributes["name"]
  e.attributes["name"] = x[0,x.length-14];
end
#------------------------------------------------------------------------------
# 3. change the material fraction references
#------------------------------------------------------------------------------
# doc.root.elements.each("//material/fraction") do |e|
#   x= e.attributes["ref"]
#   e.attributes["ref"] = x[0,x.length-14];
# end

#------------------------------------------------------------------------------
# 4. solids: names, intersections, subtractions
#------------------------------------------------------------------------------
doc.root.elements.each("//solids/*") do |e|
  x= e.attributes["name"]
  e.attributes["name"] = x[0,x.length-14];
end

doc.root.elements.each("//solids/intersection/*") do |e|
  x= e.attributes["ref"]
  if (x != nil) then
    e.attributes["ref"] = x[0,x.length-14]
  end 

  x= e.attributes["name"]
  if (x  != nil) then
    l1 = x.index('0x')
    e.attributes["name"] = x[0,l1]+x[l1+14,x.length-l1-14];
  end
end

doc.root.elements.each("//solids/subtraction/*") do |e|
  x= e.attributes["ref"]
  if (x != nil) then
    e.attributes["ref"] = x[0,x.length-14]
  end 

  x= e.attributes["name"]
  if (x  != nil) then
    l1 = x.index('0x')
    e.attributes["name"] = x[0,l1]+x[l1+14,x.length-l1-14];
  end
end

#------------------------------------------------------------------------------
# 4. structure: names
#------------------------------------------------------------------------------
doc.root.elements.each("//structure/volume") do |e|
  x= e.attributes["name"]
  e.attributes["name"] = x[0,x.length-14];

  e.elements.each("*") do |e1|
    x= e1.attributes["ref"]
    if (x != nil) then
      e1.attributes["ref"] = x[0,x.length-14]
    end
    
    x= e1.attributes["name"]
    if (x  != nil) then
      l1 = x.index('0x')
      e1.attributes["name"] = x[0,l1]+x[l1+14,x.length-l1-14];
    end
  end
end

doc.root.elements.each("//structure/volume/physvol/*") do |e|
  x= e.attributes["name"]
  if (x != nil) then
    l1 = x.index('0x')
    if (l1 >= 0) then 
      e.attributes["name"] = x[0,l1]+x[l1+14,x.length-l1-14];
    end
  end

  x= e.attributes["ref"]
  if (x != nil) then
    l1 = x.index('0x')
    if (l1 >= 0) then 
      e.attributes["ref"] = x[0,x.length-14]
    end
  end
    
end

#------------------------------------------------------------------------------
# 4. setup: ref
#------------------------------------------------------------------------------
doc.root.elements.each("//setup/*") do |e|
  x= e.attributes["ref"]
  if (x != nil) then
    l1 = x.index('0x')
    if (l1 >= 0) then 
      e.attributes["ref"] = x[0,x.length-14]
    end
  end
end

#puts root.elements[2];
doc.write(File.open("test2.gdml","w"))
