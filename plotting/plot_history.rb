# plot_profile.rb

require 'Tioga/FigureMaker'
require './plot_styles'

class ProfilePlots
  include Math
  include Tioga
  include FigureConstants
  include MyPlotStyles

  def t
    @figure_maker
  end

  def initialize
    @figure_maker = FigureMaker.default
  
    t.save_dir = 'history_plots'
    t.def_eval_function { |str| eval(str) }
  
    #@data_filename = '../LOGS-noUrca-lowT/history.data'
    @data_filename = 'history.data'
    @data_start_line = 12
    @header_array = Dvector.new(4)
    @data_array = [
      @model=Dvector.new,
      @time=Dvector.new,
      @lgMd=Dvector.new,
      @lgTeff=Dvector.new,
      @lgLsurf=Dvector.new,
      @lgLnu=Dvector.new,
      @lgLnuc=Dvector.new
    ]
    
    # observations
    @outburst_end = 54320
    @mjd = Dvector[54324,54332,54338,54371,54496,54547,54620,54752,54861]
    @obs_time = Dvector[4.36,12.21,17.65,50.90,175.73,227.12,299.71,432.47,541.49]

    @obs_Teff = Dvector[160.50,155.40,155.28,149.27,127.95,156.63,133.79,124.44,123.21]
    @derr = Dvector[3.64,2.73,1.36,1.35,5.90,2.11,2.11,3.85,1.63]
    @uerr = Dvector[3.32,2.50,1.30,1.28,4.41,1.97,1.92,3.10,1.52]
    
    @have_data = false
    @margin = 0.05
  
    t.def_enter_page_function { enter_page }
    t.def_figure('lightcurve') { lightcurve }
    t.def_figure('XTEJ1701') { compare_PR }
    t.def_figure('long-Q50') { compare_long_Q }
  end

  def enter_page
    set_default_plot_style
    # sans_serif_style
    t.default_enter_page_function
  end

  def set_boundaries(xs, ys, margin, ymin=nil, ymax=nil)
    xmin = xs.min
    xmax = xs.max
    ymin = ys.min if ymin == nil
    ymax = ys.max if ymax == nil
    width = (xmax == xmin)? 1 : xmax - xmin
    height = (ymax == ymin)? 1 : ymax - ymin
    left_boundary = xmin-margin*width
    right_boundary = xmax+margin*width
    top_boundary = ymax+margin*height
    bottom_boundary = ymin-margin*height
    return [ left_boundary, right_boundary, top_boundary, bottom_boundary ]
  end

  def read_data
    if (@have_data)
      return
    end
    Dvector.read_row(@data_filename,row=6,@header_array)
    Dvector.read(@data_filename, @data_array,start=@data_start_line)
    @have_data = true
    @grav, @Mcore, @Rcore, @Tcore = @header_array
  end

  def extract_lightcurve
    start = @lgMd.where_eq(0.0)
    kbev = 8.6173e-5
    ePhi = 0.77
    secday = 86400.0
    ys = @lgTeff[start..-1].exp10.mul(kbev*ePhi)
    xs = @time[start..-1].safe_log10(1.0)-log10(secday)
    return xs, ys
  end

  def lightcurve
    read_data
    t.do_box_labels(nil,'$\log(t/\unitstyle{d})$','$\Teff^\infty/\eV$')
    td, lt = extract_lightcurve
    xs = td.to_a
    ys = lt.to_a
    t.xaxis_log_values=true
    t.show_plot(set_boundaries([0.0,4.0],ys,@margin,ymin=0,ymax=180)) do
      t.show_polyline(xs,ys,Black,nil,LINE_TYPE_SOLID)
    end
  end
  
  def compare_PR
    t.do_box_labels(nil,'$\log(t/\unitstyle{d})$','$\Teff^\infty/\eV$')
    t.xaxis_log_values = true
    datadirs = [ '../LOGS-noUrca', '../LOGS-noUrca-lowT', '../LOGS-Urca', '../LOGS-Urca-lowT' ]
    clrs = [Black, Black, SlateGray, SlateGray ]
    stls = [LINE_TYPE_SOLID, LINE_TYPE_DASH, LINE_TYPE_SOLID, LINE_TYPE_DASH]
    thks = [ 3.0, 3.0, 1.5, 1.5 ]
    xs_arry = []
    ys_arry = []
    datadirs.size.times  do |i|
      @data_filename = datadirs[i]+'/history.data'
      @have_data = false
      read_data
      td, lt = extract_lightcurve
      xs_arry << td 
      ys_arry << lt 
    end

    xs_obs = @obs_time.log10
    ys_obs = @obs_Teff
    ue = @uerr
    de = @derr
    
    t.show_plot(set_boundaries([0.0,4.0],[0.0,180.0],@margin)) do
      datadirs.size.times do |i|
        t.line_width = thks[i]
        t.show_polyline(xs_arry[i],ys_arry[i],clrs[i],nil,stls[i])
      end
      t.show_marker('xs'=>xs_obs,'ys'=>ys_obs,'marker'=>Bullet,'scale'=>0.7, 'color'=>Blue)
      xs_obs.each_index { |i| 
        t.show_error_bars('x'=>xs_obs[i],'y'=>ys_obs[i],'dx'=>0,
        'dy_plus'=>ue[i], 'dy_minus'=>de[i], 'color'=>Blue) }
    end
  end
  
  def compare_long_Q
    t.do_box_labels(nil,'$\log(t/\unitstyle{d})$','$\Teff^\infty/\eV$')
    t.xaxis_log_values = true
    datadirs = [ '../LOGS-long-Q50', '../LOGS-long-Urca-Q50' ]
    clrs = [Black, SlateGray ]
    stls = [LINE_TYPE_SOLID, LINE_TYPE_SOLID ]
    thks = [ 3.0, 1.5 ]
    xs_arry = []
    ys_arry = []
    datadirs.size.times  do |i|
      @data_filename = datadirs[i]+'/history.data'
      @have_data = false
      read_data
      td, lt = extract_lightcurve
      xs_arry << td 
      ys_arry << lt 
    end

    xs_obs = @obs_time.log10
    ys_obs = @obs_Teff
    ue = @uerr
    de = @derr
    
    t.show_plot(set_boundaries([0.0,4.0],[100.0,200.0],@margin)) do
      datadirs.size.times do |i|
        t.line_width = thks[i]
        t.show_polyline(xs_arry[i],ys_arry[i],clrs[i],nil,stls[i])
      end
      # t.show_marker('xs'=>xs_obs,'ys'=>ys_obs,'marker'=>Bullet,'scale'=>0.7, 'color'=>Blue)
      # xs_obs.each_index { |i|
      #   t.show_error_bars('x'=>xs_obs[i],'y'=>ys_obs[i],'dx'=>0,
      #   'dy_plus'=>ue[i], 'dy_minus'=>de[i], 'color'=>Blue) }
    end
  end
  
end

ProfilePlots.new
