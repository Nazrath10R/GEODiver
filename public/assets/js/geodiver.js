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
(function() {

  GD.loadGeoDbValidation = function() {
    'use strict';
    $('#load_geo_db').validate({
      rules: {
        geo_db: { 
          geoDb: true,
          required: true
        },
      },
      submitHandler: function(form) {
        GD.loadGeoDbAjax();
      }
    });
  };



  GD.analyseValidation = function() {
    'use strict';
    $('#analyse').validate({
      rules: {},
      submitHandler: function(form) {
        GD.GeoAnalysisAjax();
      }
    });
  };

  GD.loadGeoDbAjax = function() {
    'use strict';
    var geo_db = $('input[name=geo_db]').val()
    $('#model_header_text').text('Loading GEO Dataset: ' + geo_db)
    $('#model_text').text('This should take a few seconds. Please leave this page open')
    $('#loading_modal').openModal({ dismissible: false });
    $.ajax({
      type: 'POST',
      url: '/load_geo_db',
      data: $('#load_geo_db').serialize(),
      success: function(response) {
        $('.load_geo_db_btn').removeClass('btn-large')
        $('#geo_db_summary').html(response);
        $('#geo_db_summary').show();
        $('.adv_param_collapsible').collapsible();
        $("input:radio[name=factor]:first").attr('checked', true);
        $('#' + $("input:radio[name=factor]:first").attr('id') + '_select').show()
        GD.addFactorToggle()
        $('select').material_select();
        GD.addDataSetInfo();
        GD.analyseValidation();
        $('#loading_modal').closeModal();
      },
      error: function(e, status) {
        GD.ajaxError();
      }
    });
  };

  GD.GeoAnalysisAjax = function() {
    'use strict';
    var geo_db = $('input[name=geo_db]').val()
    $('#model_header_text').text('Analysing GEO Dataset: ' + geo_db)
    $('#model_text').text('This should take a few minutes. Please leave this page open')
    $('#loading_modal').openModal({ dismissible: false});
    $.ajax({
      type: 'POST',
      url: '/analyse',
      data: $('#analyse').serialize(),
      success: function(response) {
        $('#results_section').html(response)
        $('#results_section').show()
        $('#results_tabs').tabs(); // init material tabs
        GD.createPlots()
        $('.materialboxed').materialbox(); // init materialbox
        $('#loading_modal').closeModal();
      },
      error: function(e, status) {
        GD.ajaxError();
      }
    });
  };

  GD.ajaxError = function() {
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

  GD.setUpValidatorDefaults = function() {
    $.validator.addMethod('geoDb', function (value) { 
        return /^GDS\d\d\d\d?$/.test(value); 
    }, 'Please enter a valid GEO dataset accession number (in the format GDSxxxx).');


    $.validator.setDefaults({
        errorClass: 'invalid',
        validClass: "valid",
        errorPlacement: function (error, element) {
            $(element)
                .closest("form")
                .find("label[for='" + element.attr("id") + "']")
                .attr('data-error', error.text());
        },
    });
  }

  GD.createPlots = function() {
    jsonFile = $('#DGEA').data("results-json")
    $.getJSON(jsonFile, function(json) {
      GD.createPCAPLOT(json.pc.cumVar, json.pc.expVar, json.pc.pcnames)
      GD.createVolcanoPlot(json.vol.logFC, json.vol.pVal)
    });
  }


  GD.createPCAPLOT = function(cumVar, expVar, pcaNames) {
    var CumulativePCA = { x: pcaNames, y: cumVar, type: 'scatter', name: 'Cumulative PCA' };
    var PCA = { x: pcaNames, y: expVar, type: 'scatter', name: 'PCA' };
    var data = [CumulativePCA, PCA];
    var layout = { legend: { x: 0, y: 100, traceorder: 'normal' } };
    Plotly.newPlot('PCA_plot', data, layout);
  }

  // 
  GD.createVolcanoPlot = function(xValues, yValues) {
    var trace1 = {
      x: xValues, y: yValues, mode: 'markers', type: 'scatter', name: 'volcano_plot',
      text: ['DDR1', 'RFC2', 'HSPA6', 'PAX8', 'GUCA1A', 'UBA7', 'THRA', 'PTPN21', 'CCL5', 'CYP2E1', 'EPHB3', 'ESRRA', 'CYP2A6', 'SCARB1', 'TTLL12', 'LINC00152', 'WFDC2', 'MAPK1', 'MAPK1', 'ADAM32', 'SPATA17', 'PRR22', 'PRR22', 'PXK', 'PXK', 'VPS18', 'MSANTD3', 'SLC46A1', 'SLC46A1', 'TIMD4', 'SLC39A5', 'ZDHHC11', 'ATP6V1E2', 'AFG3L1P', 'CILP2', 'CILP2', 'PIGX', 'TMEM196', 'SLC39A13', 'BEST4', 'AK9', 'CORO6', 'TMEM106A', 'TMEM106A', 'ALG10', 'ALG10', 'TTC39C', 'NEXN', 'C15orf40', 'RAX2', 'MFAP3', 'EYA3', 'GIMAP1', 'GIMAP1', 'GIMAP1', 'KLK8', 'CCDC65', 'CCDC65', 'FAM122C', 'FAM122C', 'CFAP53', 'CFAP53', 'ARMCX4', 'RBBP6', 'CENPBD1', 'TRIOBP', 'TRIOBP', 'CATSPER1', 'HOXD4', 'GSC', 'SP7', 'PDE7A', 'CNOT7', 'CRYZL1', 'PRSS33', 'PRSS33', 'C19orf26', 'C19orf26', 'MCMDC2', 'TIRAP', 'LEAP2', 'MSI2', 'SCIN', 'SCIN', 'CTCFL', 'C4orf33', 'C4orf33', 'C4orf33', 'ZNF333', 'TVP23C', 'RDH10', 'RDH10', 'SRSF12', 'FAM71A', 'FAM71A', 'GAPT', 'FLJ30901', 'ERICH5', 'ERICH5', 'CCDC185'],
      marker: { size: 7.5 }
    };
    var data = [trace1];
    var layout = {
      xaxis: { range: [-1.15, 1.15]},
      yaxis: { range: [5, 28] },
      hovermode: 'closest',
      xaxis: { title: 'Log 2 Fold Change' },
      yaxis: { title: '-Log10(P Value)' }
    };
    Plotly.newPlot('volcano_plot', data, layout);
  }

  // 
  GD.addFactorToggle = function() {
    $("input:radio[name=factor]").click(function() {
      var target = '#' + $(this).attr('id') + '_select'
      if ('#' + $('.select_factors:visible').attr('id') !== target) {
        $('.select_factors').hide()
        $(target).show()
      }
    });
  }

  GD.addDataSetInfo = function() {
    var geo_accession = $('input[name=geo_db]').val()
    var jsonFile = 'GeoDiver/DBs/' + geo_accession + '.json'
    $.getJSON(jsonFile, function(json) {
      $('#dataset_accession').text(json.Accession)
      $('#dataset_title').text(json.Title)
      $('#dataset_summary').text(json.Description)
      $('#dataset_organism').text(json.Sample_Organism)
      $('#dataset_summary').text(json.Description)
      $('#dataset_citation').text(json.Reference)
    });
  }

  GD.addUserDropDown = function() {
    $('.dropdown-button').dropdown({
      inDuration: 300, outDuration: 225, hover: true, belowOrigin: true, alignment: 'right'
    });
  }
}());




(function($) {
  $(function() {
    $('.button-collapse').sideNav();
    $('.parallax').parallax();
    $('select').material_select();
    GD.setUpValidatorDefaults()
    GD.loadGeoDbValidation();
    GD.addUserDropDown();
  });
})(jQuery);
