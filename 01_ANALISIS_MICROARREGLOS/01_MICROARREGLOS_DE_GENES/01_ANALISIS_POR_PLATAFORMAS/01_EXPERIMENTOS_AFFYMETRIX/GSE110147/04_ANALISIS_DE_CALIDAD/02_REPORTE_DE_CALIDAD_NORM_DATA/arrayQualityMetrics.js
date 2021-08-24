// (C) Wolfgang Huber 2010-2011

// Script parameters - these are set up by R in the function 'writeReport' when copying the 
//   template for this script from arrayQualityMetrics/inst/scripts into the report.

var highlightInitial = [ false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, false, false, true ];
var arrayMetadata    = [ [ "1", "IPF.1", "GSM2978752_5.CEL", "IPF", "1" ], [ "2", "IPF.2", "GSM2978753_13.CEL", "IPF", "2" ], [ "3", "IPF.3", "GSM2978754_23.CEL", "IPF", "3" ], [ "4", "IPF.4", "GSM2978755_25.CEL", "IPF", "4" ], [ "5", "IPF.5", "GSM2978756_41.CEL", "IPF", "5" ], [ "6", "IPF.6", "GSM2978757_43.CEL", "IPF", "6" ], [ "7", "IPF.7", "GSM2978759_47.CEL", "IPF", "7" ], [ "8", "IPF.8", "GSM2978760_51.CEL", "IPF", "8" ], [ "9", "IPF.9", "GSM2978761_55.CEL", "IPF", "9" ], [ "10", "IPF.10", "GSM2978762_60.CEL", "IPF", "10" ], [ "11", "IPF.11", "GSM2978764_67.CEL", "IPF", "11" ], [ "12", "IPF.12", "GSM2978765_71.CEL", "IPF", "12" ], [ "13", "IPF.13", "GSM2978766_72.CEL", "IPF", "13" ], [ "14", "IPF.14", "GSM2978767_75.CEL", "IPF", "14" ], [ "15", "IPF.15", "GSM2978768_77.CEL", "IPF", "15" ], [ "16", "IPF.16", "GSM2978770_88.CEL", "IPF", "16" ], [ "17", "IPF.17", "GSM2978771_111.CEL", "IPF", "17" ], [ "18", "IPF.18", "GSM2978772_112.CEL", "IPF", "18" ], [ "19", "IPF.19", "GSM2978773_115.CEL", "IPF", "19" ], [ "20", "Control.1", "GSM2978789_T14.CEL", "Control", "1" ], [ "21", "Control.2", "GSM2978790_T15.CEL", "Control", "2" ], [ "22", "Control.3", "GSM2978791_T16.CEL", "Control", "3" ], [ "23", "Control.4", "GSM2978792_T17.CEL", "Control", "4" ], [ "24", "Control.5", "GSM2978793_T18.CEL", "Control", "5" ], [ "25", "Control.6", "GSM2978794_T19.CEL", "Control", "6" ], [ "26", "Control.7", "GSM2978795_T20.CEL", "Control", "7" ], [ "27", "Control.8", "GSM2978796_T21.CEL", "Control", "8" ], [ "28", "Control.9", "GSM2978797_T22.CEL", "Control", "9" ], [ "29", "Control.10", "GSM2978798_T23.CEL", "Control", "10" ], [ "30", "Control.11", "GSM2978799_T24.CEL", "Control", "11" ] ];
var svgObjectNames   = [ "pca", "dens" ];

var cssText = ["stroke-width:1; stroke-opacity:0.4",
               "stroke-width:3; stroke-opacity:1" ];

// Global variables - these are set up below by 'reportinit'
var tables;             // array of all the associated ('tooltips') tables on the page
var checkboxes;         // the checkboxes
var ssrules;


function reportinit() 
{
 
    var a, i, status;

    /*--------find checkboxes and set them to start values------*/
    checkboxes = document.getElementsByName("ReportObjectCheckBoxes");
    if(checkboxes.length != highlightInitial.length)
	throw new Error("checkboxes.length=" + checkboxes.length + "  !=  "
                        + " highlightInitial.length="+ highlightInitial.length);
    
    /*--------find associated tables and cache their locations------*/
    tables = new Array(svgObjectNames.length);
    for(i=0; i<tables.length; i++) 
    {
        tables[i] = safeGetElementById("Tab:"+svgObjectNames[i]);
    }

    /*------- style sheet rules ---------*/
    var ss = document.styleSheets[0];
    ssrules = ss.cssRules ? ss.cssRules : ss.rules; 

    /*------- checkboxes[a] is (expected to be) of class HTMLInputElement ---*/
    for(a=0; a<checkboxes.length; a++)
    {
	checkboxes[a].checked = highlightInitial[a];
        status = checkboxes[a].checked; 
        setReportObj(a+1, status, false);
    }

}


function safeGetElementById(id)
{
    res = document.getElementById(id);
    if(res == null)
        throw new Error("Id '"+ id + "' not found.");
    return(res)
}

/*------------------------------------------------------------
   Highlighting of Report Objects 
 ---------------------------------------------------------------*/
function setReportObj(reportObjId, status, doTable)
{
    var i, j, plotObjIds, selector;

    if(doTable) {
	for(i=0; i<svgObjectNames.length; i++) {
	    showTipTable(i, reportObjId);
	} 
    }

    /* This works in Chrome 10, ssrules will be null; we use getElementsByClassName and loop over them */
    if(ssrules == null) {
	elements = document.getElementsByClassName("aqm" + reportObjId); 
	for(i=0; i<elements.length; i++) {
	    elements[i].style.cssText = cssText[0+status];
	}
    } else {
    /* This works in Firefox 4 */
    for(i=0; i<ssrules.length; i++) {
        if (ssrules[i].selectorText == (".aqm" + reportObjId)) {
		ssrules[i].style.cssText = cssText[0+status];
		break;
	    }
	}
    }

}

/*------------------------------------------------------------
   Display of the Metadata Table
  ------------------------------------------------------------*/
function showTipTable(tableIndex, reportObjId)
{
    var rows = tables[tableIndex].rows;
    var a = reportObjId - 1;

    if(rows.length != arrayMetadata[a].length)
	throw new Error("rows.length=" + rows.length+"  !=  arrayMetadata[array].length=" + arrayMetadata[a].length);

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = arrayMetadata[a][i];
}

function hideTipTable(tableIndex)
{
    var rows = tables[tableIndex].rows;

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = "";
}


/*------------------------------------------------------------
  From module 'name' (e.g. 'density'), find numeric index in the 
  'svgObjectNames' array.
  ------------------------------------------------------------*/
function getIndexFromName(name) 
{
    var i;
    for(i=0; i<svgObjectNames.length; i++)
        if(svgObjectNames[i] == name)
	    return i;

    throw new Error("Did not find '" + name + "'.");
}


/*------------------------------------------------------------
  SVG plot object callbacks
  ------------------------------------------------------------*/
function plotObjRespond(what, reportObjId, name)
{

    var a, i, status;

    switch(what) {
    case "show":
	i = getIndexFromName(name);
	showTipTable(i, reportObjId);
	break;
    case "hide":
	i = getIndexFromName(name);
	hideTipTable(i);
	break;
    case "click":
        a = reportObjId - 1;
	status = !checkboxes[a].checked;
	checkboxes[a].checked = status;
	setReportObj(reportObjId, status, true);
	break;
    default:
	throw new Error("Invalid 'what': "+what)
    }
}

/*------------------------------------------------------------
  checkboxes 'onchange' event
------------------------------------------------------------*/
function checkboxEvent(reportObjId)
{
    var a = reportObjId - 1;
    var status = checkboxes[a].checked;
    setReportObj(reportObjId, status, true);
}


/*------------------------------------------------------------
  toggle visibility
------------------------------------------------------------*/
function toggle(id){
  var head = safeGetElementById(id + "-h");
  var body = safeGetElementById(id + "-b");
  var hdtxt = head.innerHTML;
  var dsp;
  switch(body.style.display){
    case 'none':
      dsp = 'block';
      hdtxt = '-' + hdtxt.substr(1);
      break;
    case 'block':
      dsp = 'none';
      hdtxt = '+' + hdtxt.substr(1);
      break;
  }  
  body.style.display = dsp;
  head.innerHTML = hdtxt;
}
