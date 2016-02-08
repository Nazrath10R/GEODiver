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
        var geo_db = $('input[name=geo_db]').val()
        $('#model_header_text').text('Loading GEO Dataset: ' + geo_db)
        $('#model_text').text('This should take a few seconds. Please leave this page open')
        $('#loading_modal').openModal({ dismissible: false });
        $.ajax({
          type: 'POST',
          url: '/load_geo_db',
          data: $('#load_geo_db').serialize(),
          success: function(response) {
            $('.card-action').remove()
            $('#results_section').empty();
            $( response ).insertAfter( "#load_geo_card" );
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
      }
    });
  };

  GD.analyseValidation = function() {
    'use strict';
    $('#analyse').validate({
      rules: {},
      submitHandler: function(form) {
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
            $('.collapsible').collapsible();
            GD.createPlots()
            $('.materialboxed').materialbox(); // init materialbox
            $('#loading_modal').closeModal();
          },
          error: function(e, status) {
            GD.ajaxError();
          }
        });
      }
    });
  };

  GD.loadPcRedrawValidation = function() {
    $('#pca_redraw').validate({
      rules: {},
      submitHandler: function(form, event) {
        event.preventDefault()
        $('#principle_plot').empty()
        var x = $('select[name=PCoption1]').val()
        var y = $('select[name=PCoption2]').val()
        var jsonFile = $('#DGEA').data("results-json")
        $.getJSON(jsonFile, function(json) {
          pcaScatterPlot = GD.createPCAScatterPlot(json.pcdata, x, y)
        })
      }
    });
  }

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

  GD.createPlots = function() {
    jsonFile = $('#DGEA').data("results-json")
    $.getJSON(jsonFile, function(json) {
      pcaPlot = GD.createPCAPLOT(json.pc.cumVar, json.pc.expVar, json.pc.pcnames)
      volcanoPlot = GD.createVolcanoPlot(json.vol.logFC, json.vol.pVal, json.vol.genes)
      pcaScatterPlot = GD.createPCAScatterPlot(json.pcdata, 'PC1', 'PC2')
      GD.initialiatizePcaScatterPlot(json.pc.pcnames)
      GD.initialize_toptable(json.tops)
      $('select').material_select();
    });
    window.onresize = function() {
      Plotly.Plots.resize(pcaPlot);
      Plotly.Plots.resize(volcanoPlot);
      Plotly.Plots.resize(pcaScatterPlot);
    }
  }

  GD.createPCAScatterPlot = function(pcdata, x, y) {
    var group1 = { x: pcdata[ x + '.Group1'], y: pcdata[y + '.Group1'], type: 'scatter',mode: 'markers', name: 'Group1' };
    var group2 = { x: pcdata[ x + '.Group2'], y: pcdata[y + '.Group2'], type: 'scatter',mode: 'markers', name: 'Group2' };
    var data = [group1, group2];
    var layout = { xaxis: { title: x }, yaxis: { title: y } };

    var WIDTH_IN_PERCENT_OF_PARENT = 100
    var PCAplotGd3 = Plotly.d3.select('#principle_plot')
                       .style({width: WIDTH_IN_PERCENT_OF_PARENT + '%',
                              'margin-left': (100 - WIDTH_IN_PERCENT_OF_PARENT) / 2 + '%'});
    var pcaPlot = PCAplotGd3.node();
    Plotly.newPlot(pcaPlot, data, layout)
    return pcaPlot
  }

  GD.createPCAPLOT = function(cumVar, expVar, pcaNames) {
    var CumulativePCA = { x: pcaNames, y: cumVar, type: 'scatter', name: 'Cumulative PCA' };
    var PCA = { x: pcaNames, y: expVar, type: 'scatter', name: 'PCA' };
    var data = [CumulativePCA, PCA];
    var layout = { legend: { x: 0, y: 100, traceorder: 'normal' } };

    var WIDTH_IN_PERCENT_OF_PARENT = 100
    var PCAplotGd3 = Plotly.d3.select('#PCA_plot')
                       .style({width: WIDTH_IN_PERCENT_OF_PARENT + '%',
                              'margin-left': (100 - WIDTH_IN_PERCENT_OF_PARENT) / 2 + '%'});
    var pcaPlot = PCAplotGd3.node();
    Plotly.newPlot(pcaPlot, data, layout)
    return pcaPlot
  }

  // 
  GD.createVolcanoPlot = function(xValues, yValues, genes) {
    var trace1 = { x: xValues, y: yValues, text: genes, mode: 'markers', type: 'scatter', name: 'volcano_plot', marker: { size: 7.5 } };
    var data = [trace1];
    var layout = { xaxis: { range: [-1.15, 1.15]}, yaxis: { range: [5, 28] }, hovermode: 'closest', xaxis: { title: 'Log 2 Fold Change' }, yaxis: { title: '-Log10(P Value)' } };
    var WIDTH_IN_PERCENT_OF_PARENT = 100
    var volcanoPlotGd3 = Plotly.d3.select('#volcano_plot')
                       .style({width: WIDTH_IN_PERCENT_OF_PARENT + '%',
                              'margin-left': (100 - WIDTH_IN_PERCENT_OF_PARENT) / 2 + '%'});
    var volcanoPlot = volcanoPlotGd3.node();
    Plotly.newPlot(volcanoPlot, data, layout)
    return volcanoPlot
  }

  GD.initialize_toptable = function(dataset) {
    $('#datatable').dataTable({
      "oLanguage": {
        "sStripClasses": "",
        "sSearch": "",
        "sSearchPlaceholder": "Enter Keywords Here",
        "sInfo": "_START_ -_END_ of _TOTAL_",
        "sLengthMenu": '<span>Rows per page:</span><select class="browser-default">' +
          '<option value="10">10</option>' +
          '<option value="20">20</option>' +
          '<option value="30">30</option>' +
          '<option value="40">40</option>' +
          '<option value="50">50</option>' +
          '<option value="-1">All</option>' +
          '</select></div>'
      },
      data: dataset,
      bAutoWidth: false
    });
    $('.search-toggle').click(function() {
      if ($('.hiddensearch').css('display') == 'none')
        $('.hiddensearch').slideDown();
      else
        $('.hiddensearch').slideUp();
    });
  }

  GD.initialiatizePcaScatterPlot = function(pcnames) {
    $.each(pcnames, function(key,value) {   
      $('#PCoption1').append($("<option></option>")
        .attr("value", value).text(value));
      $('#PCoption2').append($("<option></option>")
        .attr("value", value).text(value));
    });
    GD.loadPcRedrawValidation()
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
