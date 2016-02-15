

// GD module
(function() {
  GD.download_all_results = function () {
    $('#history_container').on('click', '#download-all-results', function(event) {
      event.preventDefault();
      $('#model_header_text').text('Creating Download Link');
      $('#model_text').text('This should take a few seconds. Please leave this page open');
      $('#loading_modal').openModal({ dismissible: false});
      $.fileDownload($(this).data('download'), {
          successCallback: function(url) {
            $('#loading_modal').closeModal();
          },
          failCallback: function(responseHtml, url) {
            $('#loading_modal').closeModal();
          }
      });
      $('#loading_modal').closeModal();
      return false; //this is critical to stop the click event which will trigger a normal file download!
    });
  };

  GD.delete_result = function () {
    $('#history_container').on('click', '#delete_results', function(event) {
      $('#model_header_text').text('Deleting Result');
      $('#model_text').text('This should take a few seconds. Please leave this page open');
      $('#loading_modal').openModal({ dismissible: false});
      var resultId = $(this).closest('.card').data('result');
      var geoDb =  $(this).closest('.card').data('geodb');
      $.ajax({
        type: 'POST',
        url: '/delete_result',
        data: {result_id: resultId, geo_db: geoDb},
        success: function(response) {
          location.reload();
        },
        error: function(e, status) {
          GD.ajaxError(e, status);
        }
      });
    });
  };

  GD.share_result = function () {
    $('#history_container').on('click', '#share_btn', function() {
      var share_link =  $(this).closest('.card').data('share-link');
      $('#share_the_link_btn').show();
      $('#share_btn').hide();
      $('#share_link_input').val(share_link);
      $('#share_link_input').prop("readonly", true);
      $('#share_model').openModal();
      $('#share_model').attr('data-share-link', share_link);
      $('#share_link_input').select();
      $.ajax({
        type: 'POST',
        url: share_link,
        error: function(e, status) {
          GD.ajaxError(e, status);
        }
      });
    });
    $('#history_container').on('click', '#share_the_link_btn', function() {
      var share_link =  $(this).closest('.card').data('share-link');
      $('#share_link_input1').val(share_link);
      $('#share_link_input1').prop("readonly", true);
      $('#share_the_link_model').openModal();
      $('#share_the_link_model').attr('data-share-link', share_link);
      $('#share_link_input1').select();
    });

    $(".share_link_input").focus(function() {
      $(this).select();
      // Work around Chrome's little problem
      $(this).mouseup(function() {
          // Prevent further mouseup intervention
          $(this).unbind("mouseup");
          return false;
      });
    });
  };

  GD.remove_share = function () {
    $('#history_container').on('click', '.remove_link', function(event) {
      var share_link =  $(this).closest('.modal').data('share-link');
      var remove_link = share_link.replace(/\/sh\//, '/rm/');
      $('#share_the_link_btn').hide();
      $('#share_btn').show();
      $('#share_model').closeModal();
      $('#share_the_link_model').closeModal();
      $.ajax({
        type: 'POST',
        url: remove_link,
        error: function(e, status) {
          GD.ajaxError(e, status);
        }
      });
    });
  };


}());



$( document ).ready(function() {
  $('.modal-trigger').leanModal();
  GD.download_all_results();
  GD.delete_result();
  GD.share_result();
  GD.remove_share();
});
