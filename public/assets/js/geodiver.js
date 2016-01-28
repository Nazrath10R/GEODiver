/*
    GD - GeoDiver's JavaScript module

    Define a global GD (acronym for GeoDiver) object containing all
    GD associated methods:
*/

//define global GD object
var GD;
if (!GD) {
    GD = {};
}

//GD module
(function () {
  GD.createPlots = function () {
    jsonFile = $('#DGEA').data("results-json")
    $.getJSON(jsonFile, function( json ) {
      GD.createPCAPLOT(json.pc.cumVar, json.pc.expVar, json.pc.pcnames)
      GD.createVolcanoPlot(json.vol.logFC, json.vol.pVal)
    });
  }


  GD.createPCAPLOT = function(cumVar, expVar, pcaNames) {
    var CumulativePCA = {
      x: pcaNames,
      y: cumVar,
      type: 'scatter',
      name: 'Cumulative PCA'
    };
    var PCA = {
      x: pcaNames, 
      y: expVar, 
      type: 'scatter',
      name: 'PCA'
    };
    var data = [CumulativePCA, PCA];
    var layout = {legend: {
      x: 0,
      y: 100,
      traceorder: 'normal',
      font: {
        family: 'sans-serif',
        size: 12,
        color: '#000'
      },
      bgcolor: '#E2E2E2',
      bordercolor: '#FFFFFF',
      borderwidth: 2
    }};
    Plotly.newPlot('PCA_plot', data, layout);
  }

  // 
  GD.createVolcanoPlot = function(xValues, yValues){
    var trace1 = {
      // x: [-1.065,-1.051,-1.004,-1.000,-0.999,-0.980,-0.973,-0.973,0.946,-0.945,-0.934,-0.933,-0.926,-0.922,-0.911,-0.910,-0.907,-0.904,-0.899,-0.888,-0.884,0.879,-0.872,-0.865,-0.865,-0.862,-0.862,-0.853,-0.848,-0.846,-0.832,0.831,-0.831,-0.830,-0.829,-0.825,-0.823,-0.817,-0.813,0.811,-0.811,-0.810,-0.806,-0.805,0.796,0.790,-0.790,-0.784,-0.783,0.778,0.775,0.768,-0.767,-0.765,-0.764,0.762,-0.762,-0.756,-0.756,-0.755,-0.753,-0.752,-0.752,0.752,0.750,-0.746,-0.738,-0.728,0.727,-0.726,0.725,-0.724,-0.723,-0.722,0.720,-0.719,0.715,-0.713,-0.712,-0.711,-0.710,-0.710,-0.709,-0.709,-0.708,0.708,-0.708,-0.705,-0.704,-0.704,0.704,-0.702,-0.701,-0.700,-0.697,-0.697,-0.697,-0.696,-0.694,-0.693],
      // y: [18.719, 25.497, 16.672, 20.371, 18.478, 16.587, 17.416, 24.482, 16.318, 14.420, 15.494, 19.096, 18.949, 19.273, 15.724, 20.303, 13.869, 16.857, 14.971, 18.154, 15.202, 13.980, 13.242, 17.531, 12.710, 14.716, 17.522, 17.252, 16.381, 16.474, 16.152, 16.055, 14.314, 14.649, 16.572, 12.504, 15.049, 18.970, 16.402, 11.053, 17.033, 12.371, 16.247, 16.483, 9.519, 13.128, 15.610, 16.049, 13.173, 12.177, 13.406, 14.160, 14.090, 18.330, 15.442, 11.627, 16.060, 20.480, 15.637, 10.734, 12.833, 14.146, 16.210, 12.845, 12.499, 16.317, 14.917, 15.956, 13.772, 13.377, 10.302, 10.077, 14.405, 15.568, 11.233, 16.008, 12.924, 12.902, 15.783, 13.669, 16.335, 12.782, 15.620, 14.087, 10.620, 14.584, 14.576, 10.031, 11.019, 16.549, 13.658, 10.919, 14.061, 7.300, 13.507, 13.492, 13.230, 16.277, 12.482, 8.150],
      x: xValues,
      y: yValues,
      mode: 'markers',
      type: 'scatter',
      name: 'volcano_plot',
      text: ['DDR1','RFC2','HSPA6','PAX8','GUCA1A','UBA7','THRA','PTPN21','CCL5','CYP2E1','EPHB3','ESRRA','CYP2A6','SCARB1','TTLL12','LINC00152','WFDC2','MAPK1','MAPK1','ADAM32','SPATA17','PRR22','PRR22','PXK','PXK','VPS18','MSANTD3','SLC46A1','SLC46A1','TIMD4','SLC39A5','ZDHHC11','ATP6V1E2','AFG3L1P','CILP2','CILP2','PIGX','TMEM196','SLC39A13','BEST4','AK9','CORO6','TMEM106A','TMEM106A','ALG10','ALG10','TTC39C','NEXN','C15orf40','RAX2','MFAP3','EYA3','GIMAP1','GIMAP1','GIMAP1','KLK8','CCDC65','CCDC65','FAM122C','FAM122C','CFAP53','CFAP53','ARMCX4','RBBP6','CENPBD1','TRIOBP','TRIOBP','CATSPER1','HOXD4','GSC','SP7','PDE7A','CNOT7','CRYZL1','PRSS33','PRSS33','C19orf26','C19orf26','MCMDC2','TIRAP','LEAP2','MSI2','SCIN','SCIN','CTCFL','C4orf33','C4orf33','C4orf33','ZNF333','TVP23C','RDH10','RDH10','SRSF12','FAM71A','FAM71A','GAPT','FLJ30901','ERICH5','ERICH5','CCDC185'],
      marker: { size: 7.5 }
    };
    var data = [ trace1 ];
    var layout = {
      xaxis: { range: [ -1.15, 1.15 ] },
      yaxis: { range: [5, 28] },
      hovermode: 'closest',
      xaxis: {
        title: 'Log 2 Fold Change',
        titlefont: { family: 'Courier New, monospace', size: 18, color: '#7f7f7f' }
      },
      yaxis: {
        title: '-Log10(P Value)',
        titlefont: { family: 'Courier New, monospace', size: 18, color: '#7f7f7f' }
      }
     };
    Plotly.newPlot('volcano_plot', data, layout);
  }

  // 
  GD.addFactorToggle = function () {
    $("input:radio[name=factor]").click(function(){
      var target = '#' + $(this).attr('id') + '_select'
      if ( '#' + $('.select_factors:visible').attr('id') !== target ) {
        $('.select_factors').hide()      
        $(target).show()
      }
    });
  }

  GD.addDataSetInfo = function () {
    var geo_accession = $('input[name=geo_db]').val()
    var jsonFile = 'GeoDiver/DBs/' + geo_accession + '.json'
    $.getJSON(jsonFile, function( json ) {
      $('#dataset_accession').text(json.Accession)
      $('#dataset_title').text(json.Title)
      $('#dataset_summary').text(json.Description)
      $('#dataset_organism').text(json.Sample_Organism)
      $('#dataset_summary').text(json.Description)
      $('#dataset_citation').text(json.Reference)
    });
  }

  GD.addUserDropDown =function () {
    $('.dropdown-button').dropdown({
      inDuration: 300,
      outDuration: 225,
      hover: true, // Activate on hover
      belowOrigin: true, // Displays dropdown below the button
      alignment: 'right' // Displays dropdown with edge aligned to the left of button
    });
  }
}());

// A function that validates the input - Utilises Jquery.Validator.js
var inputValidation = function () {
  'use strict';
  $('#input').validate({
    rules: {
        input: {  // Update id and block 
            minlength: 5,
            maxlength: 7,
            required: true,
        },
    },
    errorClass: 'invalid',
    errorPlacement: function (error, element) {
        element.next("label").attr("data-error", error.contents().text());
    },
    submitHandler: function(form) {
      var geo_db = $('input[name=geo_db]').val()
      if ($('#input').attr('action') === '/load_geo_db') {
        $('#model_header_text').text('Loading GEO Dataset: ' + geo_db)
        $('#model_text').text('This should take a few seconds. Please leave this page open')
      } else if ($('#input').attr('action') === '/') {
        $('#model_header_text').text('Analysing GEO Dataset: ' + geo_db)
        $('#model_text').text('This should take a few minutes. Please leave this page open')
      }
      $('#loading_modal').openModal({dismissible: false});
      ajaxFunction();
    }
  });
};


// Sends the data within the form to the Server
var ajaxFunction = function () {
  'use strict';
  $.ajax({
    type: 'POST',
    url: $('#input').attr('action'),
    data: $('#input').serialize(),
    success: function(response) {
      if ($('#input').attr('action') === '/load_geo_db') {
        $('#geo_db_summary').html(response);
        $('#geo_db_summary').show();
        $('.adv_param_collapsible').collapsible();
        $('#analyse_geo_btn').html('Analyse<i class="material-icons right">send</i>')
        $("input:radio[name=factor]:first").attr('checked', true);        
        $('#'+ $("input:radio[name=factor]:first").attr('id') + '_select').show()
        GD.addFactorToggle()
        $('select').material_select();
        GD.addDataSetInfo();
        $('#input').attr('action', '/analyse')
        $('#loading_modal').closeModal();
      } else if ($('#input').attr('action') === '/analyse') {
        $('#results_section').html(response)
        $('#results_section').show()
        $('.materialboxed').materialbox(); // init materialbox
        $('#results_tabs').tabs(); // init material tabs
        GD.createPlots()
        $('#loading_modal').closeModal();
      }
    },
    error: function (e, status) {
      var errorMessage;
      if (e.status == 500 || e.status == 400) {
        errorMessage = e.responseText;
        $('#results_section').show();
        $('#results_section').html(errorMessage);
        $('#loading_modal').closeModal(); // remove progress notification
      } else {
        errorMessage = e.responseText;
        $('#results_section').show();
        $('#results_section').html('There seems to be an unidentified Error.');
        $('#loading_modal').closeModal(); // remove progress notification
      }
    }
  });
};


(function($){
  $(function(){

    $('.button-collapse').sideNav();
    $('.parallax').parallax();
    $('select').material_select();
    inputValidation();
    GD.addUserDropDown();

  });
})(jQuery);
