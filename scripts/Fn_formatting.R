# functions for formatting the cancer site names

fn_format.site=function(site.names){
  site.names=gsub(' and other nervous system',"/CNS",site.names)
  site.names=gsub(' and renal pelvis',"",site.names)
  site.names=gsub('other endocrine including thymus',"other\nendocrine",site.names)
  site.names=gsub(' including heart',"",site.names)
  site.names=gsub('colon and rectum',"colorectal",site.names)
  site.names=gsub("corpus and uterus  nos","uterus",site.names)
  site.names=gsub('non-hodgkin',"NH",site.names)
  site.names=gsub(' and intrahepatic bile duct',"/IBL",site.names)
  site.names=gsub(' and pharynx',"",site.names)
  site.names=gsub(' of the skin',"",site.names)
  site.names=gsub(' including thymus',"",site.names)
  site.names=gsub(' including heart',"",site.names)
  site.names=gsub('  anal canal and anorectum',"",site.names)
  site.names=gsub('lung and bronchus',"lung",site.names)
  site.names=gsub('urinary bladder',"bladder",site.names)
  site.names=gsub('small intestine',"small\nintestine",site.names)
  site.names=gsub(' organs',"",site.names)
  
  site.names
}
fn_format.site2=function(site.names){
  # site.names=gsub(' and other nervous system',"/CNS",site.names)
  # site.names=gsub(' and renal pelvis',"",site.names)
  # site.names=gsub('other endocrine including thymus',"other\nendocrine",site.names)
  # site.names=gsub(' including heart',"",site.names)
  site.names=gsub('colon and rectum',"colorectal",site.names)
  site.names=gsub("corpus and uterus  nos","corpus and uterus, nos",site.names)
  site.names=gsub('non hodgkin lymphoma',"non-hodgkin's lymphoma",site.names)
  # site.names=gsub(' and intrahepatic bile duct',"/IBL",site.names)
  # site.names=gsub(' and pharynx',"",site.names)
  # site.names=gsub(' of the skin',"",site.names)
  # site.names=gsub(' including thymus',"",site.names)
  # site.names=gsub(' including heart',"",site.names)
  site.names=gsub('  anal canal and anorectum',", anal canal and anorectum",site.names)
  # site.names=gsub('lung and bronchus',"lung",site.names)
  site.names=gsub('urinary bladder',"bladder",site.names)
  # site.names=gsub('small intestine',"small\nintestine",site.names)
  site.names=gsub('nose ',"nose,",site.names)
  site.names=gsub('non epithelial','non-epithelial',site.names)
  site.names=gsub('trachea ',"trachea,",site.names)
  site.names=gsub('peritoneum ','peritoneum,',site.names)
  site.names=gsub(' organs',"",site.names)
  
  
  site.names
}
