/*
    GD - GeoDiver's JavaScript module

    Define a global GD (acronym for GeoDiver) object containing all
    GD associated methods:
*/

// define global GD object
var GD;
if (!GD) {
  GD = {};
}

// GD module
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
  };

  ///// AJAXs

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
        var geo_db = $('input[name=geo_db]').val();
        $('#model_header_text').text('Loading GEO Dataset: ' + geo_db);
        $('#model_text').text('This should take a few seconds. Please leave this page open');
        $('#loading_modal').openModal({ dismissible: false });
        $.ajax({
          type: 'POST',
          url: '/load_geo_db',
          data: $('#load_geo_db').serialize(),
          success: function(response) {
            $('.card-action').remove();
            $('#results_section').empty();
            $( response ).insertAfter( "#load_geo_card" );
            $('#geo_db_summary').html(response);
            $('#geo_db_summary').show();
            $('.adv_param_collapsible').collapsible();
            $("input:radio[name=factor]:first").attr('checked', true);
            $('#' + $("input:radio[name=factor]:first").attr('id') + '_select').show();
            GD.addFactorToggle();
            $('select').material_select();
            GD.addDataSetInfo();
            GD.analyseValidation();
            $('#loading_modal').closeModal();
          },
          error: function(e, status) {
            GD.ajaxError(e, status);
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
        var geo_db = $('input[name=geo_db]').val();
        $('#model_header_text').text('Analysing GEO Dataset: ' + geo_db);
        $('#model_text').text('This should take a few minutes. Please leave this page open');
        $('#loading_modal').openModal({ dismissible: false});
        $.ajax({
          type: 'POST',
          url: '/analyse',
          data: $('#analyse').serialize(),
          success: function(response) {
            $('#results_section').html(response);
            $('#results_section').show();
            $('#results_tabs').tabs(); // init material tabs
            GD.createPlots();
            $('.materialboxed').materialbox(); // init materialbox
            $('#loading_modal').closeModal();
          },
          error: function(e, status) {
            GD.ajaxError(e, status);
          }
        });
      }
    });
  };

  GD.geneExpressionAjax = function(currentRow, geneId) {
    $('#model_header_text').text('Loading Graphics for Gene: ' + geneId);
    $('#model_text').text('This should take a few seconds. Please leave this page open');
    $('#loading_modal').openModal({ dismissible: false});
    var resultId = currentRow.closest('.results_card').data('results_id');
    var geoDb = currentRow.closest('.results_card').data('geo_db');
    console.log(geoDb);
    console.log(resultId);
    console.log(geneId);
    console.log('dffdfd');

    $.ajax({
      type: 'POST',
      url: '/gene_expression_url',
      data: {gene_id: geneId, result_id: resultId, geo_db: geoDb},
      success: function(response) {
        currentRow.addClass('parent');
        currentRow.after('<tr class="child" id="'+ geneId + 'ChildRow"><td colspan="8"><div id="' + geneId + 'Plot"></div></td></tr>');
        GD.createExpressionPlot(response, geneId);
        $('#loading_modal').closeModal();
      },
      error: function(e, status) {
        GD.ajaxError(e, status);
      }
    });
  };

  GD.interactionNetworkAjax = function(currentRow, pathId) {
    $('#model_header_text').text('Loading Graphics for GeneSet: ' + pathId);
    $('#model_text').text('This should take a few seconds. Please leave this page open');
    $('#loading_modal').openModal({ dismissible: false});
    var resultId = currentRow.closest('.results_card').data('results_id');
    var geoDb = currentRow.closest('.results_card').data('geo_db');
    console.log(geoDb);
    console.log(resultId);
    console.log(pathId);

    $.ajax({
      type: 'POST',
      url: '/interaction',
      data: {path_id: pathId, result_id: resultId, geo_db: geoDb},
      success: function(response) {
        console.log(response);
        currentRow.addClass('parent');
        currentRow.after('<tr class="child" id="'+ pathId + 
                         'ChildRow"><td colspan="8"><div id="' + pathId + 'Plot">' +
                         '<p>hey</p>' +
                         '</div></td></tr>');
        $('#loading_modal').closeModal();
      },
      error: function(e, status) {
        GD.ajaxError(e, status);
      }
    });  };

  GD.ajaxError = function(e, status) {
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
  };

  GD.loadPcRedrawValidation = function() {
    $('#pca_redraw').validate({
      rules: {},
      submitHandler: function(form, event) {
        event.preventDefault(); // because we're not submitting at all
        $('#principle_plot').empty();
        var x = $('select[name=PCoption1]').val();
        var y = $('select[name=PCoption2]').val();
        var jsonFile = $('#overview').data("overview-json");
        $.getJSON(jsonFile, function(json) {
          GD.createPCAScatterPlot(json.pcdata, x, y);
        });
      }
    });
  };


  ///// PLOTS

  GD.createPlots = function() {
    var jsonFile = $('#overview').data("overview-json");
    $.getJSON(jsonFile, function(json) {
      pcaPlot = GD.createPCAPLOT(json.pc.cumVar, json.pc.expVar, json.pc.pcnames);
      pcaScatterPlot = GD.createPCAScatterPlot(json.pcdata, 'PC1', 'PC2');
      GD.initialiatizePcaScatterPlot(json.pc.pcnames);
      $('select').material_select();
    });

    var dgeaJsonFile = $('#DGEA').data("dgea-json");
    $.getJSON(dgeaJsonFile, function(json) {
      volcanoPlot = GD.createVolcanoPlot(json.vol.logFC, json.vol.pVal, json.vol.genes);
      GD.initializeToptable(json.tops, 'dgea-top-table', 'dgea-top-table-wrapper');
    });

    var gseaJsonFile = $('#GSEA').data("gsea-json");
    $.getJSON(gseaJsonFile, function(json) {
      GD.initializeToptable(json.tops, 'gsea-top-table', 'gsea-top-table-wrapper');
    });
    window.onresize = function() {
      Plotly.Plots.resize(pcaPlot);
      Plotly.Plots.resize(volcanoPlot);
      Plotly.Plots.resize(pcaScatterPlot);
    };
  };

  GD.createPCAScatterPlot = function(pcdata, x, y) {
    var group1 = { x: pcdata[ x + '.Group1'], y: pcdata[y + '.Group1'], type: 'scatter', mode: 'markers', name: 'Group1' };
    var group2 = { x: pcdata[ x + '.Group2'], y: pcdata[y + '.Group2'], type: 'scatter', mode: 'markers', name: 'Group2' };
    var data = [group1, group2];
    var layout = { xaxis: { title: x }, yaxis: { title: y } };

    var parentWidth = 100;
    var PCAplotGd3 = Plotly.d3.select('#principle_plot')
                       .style({width: parentWidth + '%',
                              'margin-left': (100 - parentWidth) / 2 + '%'});
    var pcaPlot = PCAplotGd3.node();
    Plotly.newPlot(pcaPlot, data, layout);
    return pcaPlot;
  };

  GD.createPCAPLOT = function(cumVar, expVar, pcaNames) {
    var CumulativePCA = { x: pcaNames, y: cumVar, type: 'scatter', name: 'Cumulative PCA' };
    var PCA = { x: pcaNames, y: expVar, type: 'scatter', name: 'PCA' };
    var data = [CumulativePCA, PCA];
    var layout = { legend: { x: 0, y: 100, traceorder: 'normal' } };

    var parentWidth = 100;
    var PCAplotGd3 = Plotly.d3.select('#PCA_plot')
                       .style({width: parentWidth + '%',
                              'margin-left': (100 - parentWidth) / 2 + '%'});
    var pcaPlot = PCAplotGd3.node();
    Plotly.newPlot(pcaPlot, data, layout);
    return pcaPlot;
  };

  // 
  GD.createVolcanoPlot = function(xValues, yValues, genes) {
    var trace1 = { x: xValues, y: yValues, text: genes, mode: 'markers', type: 'scatter', name: 'volcano_plot', marker: { size: 7.5 } };
    var data = [trace1];
    var layout = { xaxis: { title: 'Log 2 Fold Change'}, yaxis: { title: '-Log10(P Value)' }, hovermode: 'closest' };
    var parentWidth = 100;
    var volcanoPlotGd3 = Plotly.d3.select('#volcano_plot')
                       .style({width: parentWidth + '%',
                              'margin-left': (100 - parentWidth) / 2 + '%'});
    var volcanoPlot = volcanoPlotGd3.node();
    Plotly.newPlot(volcanoPlot, data, layout);
    return volcanoPlot;
  };

  GD.createExpressionPlot = function (response, geneId) {
    var trace1 = { x: response.group1.x, y: response.group1.y, type: 'bar', name: 'Group 1' };
    var trace2 = { x: response.group2.x, y: response.group2.y, type: 'bar', name: 'Group 2' };
    var data = [trace1, trace2];
    var layout = { barmode: 'group', xaxis: { title: 'Sample', tickangle: -40, position: -0.5}, yaxis: { title: 'Expression' } };
    var parentWidth = 100;
    var expressionPlotGd3 = Plotly.d3.select('#' + geneId + 'Plot')
                       .style({width: parentWidth + '%',
                              'margin-left': (100 - parentWidth) / 2 + '%'});
    var expressionPlot = expressionPlotGd3.node();
    Plotly.newPlot(expressionPlot, data, layout);
    window.onresize = function() {
      Plotly.Plots.resize(expressionPlot);
    };
  };


  ///// General GD functions

  GD.initializeToptable = function(dataset, tableId, tableWrapperId) {
    dataset = GD.addPlotIconToTopTable(dataset);
    var dt = $('#' + tableId).dataTable({
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
      order: [[ 5, "asc" ]],
      bAutoWidth: false
    });

    GD.makePlotIconClickable(tableWrapperId);

    $('#' + tableWrapperId).on('click', '.search-toggle', function() {
      if ($('.hiddensearch').css('display') == 'none')
        $('.hiddensearch').slideDown();
      else
        $('.hiddensearch').slideUp();
    });

    $('#' + tableWrapperId).on('click', '.download-top-table', function() {
        $('#model_header_text').text('Creating Download Link');
        $('#model_text').text('This should take a few seconds. Please leave this page open');
        $('#loading_modal').openModal({ dismissible: false});
        $.fileDownload($(this).attr('href'), {
            successCallback: function(url) {
              $('#loading_modal').closeModal();
            },
            failCallback: function(responseHtml, url) {
              $('#loading_modal').closeModal();
            }
        });
        return false; //this is critical to stop the click event which will trigger a normal file download!
    });
  };

  GD.addPlotIconToTopTable = function(dataset) {
    $.each(dataset, function(key, row) {
      row.push( '<i class="material-icons child-row-chart light-blue-text text-darken-3">insert_chart</i>' );
    });
    return (dataset);
  };

  GD.makePlotIconClickable = function(tableWrapperId) {
    $('#' + tableWrapperId).on('click', '.child-row-chart', function() {
      var currentRow = $(this).closest('tr');
      var geneId = currentRow.children('td:first').text();
      if ( $( '#' + geneId + 'ChildRow' ).length ) {
        $('#' + geneId + 'ChildRow').remove();
      } else {
        if (tableWrapperId === 'dgea-top-table-wrapper') {
          GD.geneExpressionAjax(currentRow, geneId);
        } else if (tableWrapperId === 'gsea-top-table-wrapper'){
          GD.interactionNetworkAjax(currentRow, geneId);
        }
      }
    });
  };

  GD.initialiatizePcaScatterPlot = function(pcnames) {
    $.each(pcnames, function(key,value) {   
      $('#PCoption1').append($("<option></option>")
        .attr("value", value).text(value));
      $('#PCoption2').append($("<option></option>")
        .attr("value", value).text(value));
    });
    GD.loadPcRedrawValidation();
  };

  // 
  GD.addFactorToggle = function() {
    $("input:radio[name=factor]").click(function() {
      var target = '#' + $(this).attr('id') + '_select';
      if ('#' + $('.select_factors:visible').attr('id') !== target) {
        $('.select_factors').hide();
        $(target).show();
      }
    });
  };

  GD.addDataSetInfo = function() {
    var geoAccession = $('input[name=geo_db]').val();
    var jsonFile = 'GeoDiver/DBs/' + geoAccession + '.json';
    $.getJSON(jsonFile, function(json) {
      $('#dataset_accession').text(json.Accession);
      $('#dataset_title').text(json.Title);
      $('#dataset_summary').text(json.Description);
      $('#dataset_organism').text(json.Sample_Organism);
      $('#dataset_summary').text(json.Description);
      $('#dataset_citation').text(json.Reference);
    });
  };

  GD.addUserDropDown = function() {
    $('.dropdown-button').dropdown({
      inDuration: 300, outDuration: 225, hover: true, belowOrigin: true, alignment: 'right'
    });
  };
}());


(function($) {
  $(function() {
    $('.button-collapse').sideNav();
    $('.parallax').parallax();
    $('select').material_select();
    GD.setUpValidatorDefaults();
    GD.loadGeoDbValidation();
    GD.addUserDropDown();
  });
})(jQuery);
