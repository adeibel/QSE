# plot_profile.rb

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

    t.save_dir = 'profile_plots'
    t.def_eval_function { |str| eval(str) }

   # @data_dir = '../LOGS-long-Urca-Q50'
   	@data_dir = ''
    @profile_basename = 'profile'
    @profile_id = 47
    @data_start_line = 12
    @header_array = Dvector.new(7)
    @data_array = [
      @zone=Dvector.new,
      @mass=Dvector.new,
      @dm=Dvector.new,
      @area=Dvector.new,
      @gravity=Dvector.new,
      @eLambda=Dvector.new,
      @ePhi=Dvector.new,
      @temperature=Dvector.new,
      @luminosity=Dvector.new,
      @pressure=Dvector.new,
      @density=Dvector.new,
      @zbar=Dvector.new,
      @abar=Dvector.new,
      @Xn=Dvector.new,
      @Qimp=Dvector.new,
      @cp=Dvector.new,
      @eps_nu=Dvector.new,
      @eps_nuc=Dvector.new,
      @kcond=Dvector.new
    ]

    @history_header_array = Dvector.new(4)
    @history_array = [
      @model=Dvector.new,
      @time=Dvector.new,
      @lgMd=Dvector.new,
      @lgTeff=Dvector.new,
      @lgLsurf=Dvector.new,
      @lgLnu=Dvector.new,
      @lgLnuc=Dvector.new
    ]

    @have_data = false
    @margin = 0.05

    t.def_enter_page_function { enter_page }
    
    t.def_figure('temperature') { temperature }
    t.def_figure('heating_profiles') { temperature_profiles }
    t.def_figure('cooling_profiles') { cooling_profiles }
    t.def_figure('thermal_timescale') { thermal_timescale(@profile_id) }
    t.def_figure('neutrino_timescale') { neutrino_timescale(@profile_id) }
    t.def_figure('both_timescales') { show_both_timescales_with_legend(@profile_id) }
    # t.def_figure('density') { density }
    # t.def_figure('specific_heat') { specific_heat }
    # t.def_figure('conductivity') { conductivity }
    # t.def_figure('composition') { comp_with_zones }
    # t.def_figure('enu') { enu_with_zones }
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

  def read_data(data_filename)
    Dvector.read_row(data_filename,row=6,@header_array)
    Dvector.read(data_filename, @data_array,start=@data_start_line)
  end

  def read_history(datadir)
    datafile = datadir + '/history.data'
    Dvector.read(datafile, @history_array,start=@data_start_line)
  end

  def parse_history(datadir,heating=true)
    read_history(datadir)
    # trim to start where t ~ 1 hr and then go logarithmically up to end of outburst
    lgtd = @time.div(86400.0).safe_log10
    outburst_end = @lgMd.where_last_gt(0.0)
    if heating then
      first_point = lgtd.where_first_ge(-1.0)
      last_point = outburst_end
    else
      first_point = outburst_end+lgtd[outburst_end+1..-1].where_first_ge(-1.0)
      # print "first point = ",first_point,"\n"
      last_point = -1
    end
    profiles = []
    current_lgtd = lgtd[first_point]
    m_last = @model[last_point]
    lgtd[first_point..last_point].each2(@model[first_point..last_point]) do |p,m|
      if p > current_lgtd + 0.5 or m == m_last then
        profiles << m.to_int
        current_lgtd = p
      end
    end
    return profiles
  end

  def temperature_profiles
    t.do_box_labels(nil,'$P/g$ [\columnunit]','$T$ [$\Giga\K$]')
    t.xaxis_log_values = true
    profile_ids = parse_history(@data_dir)
    xs_arry = []
    ys_arry = []
    time = []
    profile_ids.each do |p|
      data_filename = @data_dir + '/profile'+sprintf('%0.4d',p)
      read_data(data_filename)
      time << @header_array[1].div(86400.0)
      ys_arry << @temperature.mul(1.0e-9)
      xs_arry << @pressure.div(@gravity).log10
      # print "time = "+sprintf("%11.4f",time[-1])+"\n"
    end
    y_lab = 0.15
    t.show_plot(set_boundaries(xs_arry[0], ys_arry[0], @margin, ymin = ys_arry[0].min, ymax = ys_arry[-1].max)) do
      xs_arry.size.times do |i|
        lc = (i < xs_arry.size-1) ? SlateGray : Black 
        t.show_polyline(xs_arry[i],ys_arry[i],lc,nil,LINE_TYPE_SOLID)
        x_lab = xs_arry[i][ys_arry[i].where_closest(y_lab)]
        t.show_label('at'=>[x_lab,y_lab],'text'=>sprintf('%d d',time[i]),'color'=>lc,'justification'=>LEFT_JUSTIFIED)
        y_lab += 0.05
      end
    end
  end

  def cooling_profiles
    t.do_box_labels(nil,'$P/g$ [\columnunit]','$T$ [$\Giga\K$]')
    t.xaxis_log_values = true
    profile_ids = parse_history(@data_dir,heating=false)
    xs_arry = []
    ys_arry = []
    time = []
    profile_ids.each do |p|
      data_filename = @data_dir + '/profile'+sprintf('%0.4d',p)
      read_data(data_filename)
      time << @header_array[1].div(86400.0)
      ys_arry << @temperature.mul(1.0e-9)
      xs_arry << @pressure.div(@gravity).log10
    end
    xlim = [xs_arry[0].min-0.75,xs_arry[0].max]
    t.show_plot(set_boundaries(xlim, ys_arry[0], @margin, ymin = ys_arry[-1].min, ymax = ys_arry[0].max)) do
      xs_arry.size.times do |i|
        lc = (i < xs_arry.size-1) ? SlateGray : Black 
        t.show_polyline(xs_arry[i],ys_arry[i],lc,nil,LINE_TYPE_SOLID)
        x_lab = xs_arry[i][0]
        y_lab = ys_arry[i][0]
        t.show_label('at'=>[x_lab,y_lab],'text'=>sprintf('%d d ',time[i]),'color'=>lc,'justification'=>RIGHT_JUSTIFIED)
        y_lab -= 0.05
      end
    end
  end

  def temperature(profile_id)
    data_filename = @data_dir + '/profile'+sprintf('%0.4d',profile_id)
    print '\n\n'
    print data_filename
    read_data(data_filename)
    t.do_box_labels(nil,'$P/g$ [\columnunit]','$T$ [$\Giga\K$]')
    ys = @temperature.mul(1.0e-9)
    xs = @pressure.div(@gravity).log10
    t.xaxis_log_values = true
    t.show_plot(set_boundaries(xs,ys,@margin)) do
      t.show_polyline(xs,ys)
      # t.show_text('at'=>[13.5,9.5],'text'=>sprintf('%11.4e',@time)+' s')
    end
  end

  def eval_thermal_timescale
    chi = @cp.div(@density).div(@kcond).sqrt
    dtau = chi.mul(@dm).div(@area)
    tau = Dvector.new(@dm.size)
    tau[0] = sqrt(@cp[0]/@density[0]/@kcond[0])*(@pressure[0]/@gravity[0])
    dtau.size.times do |i|
      tau[i] = tau[i-1] + dtau[i] unless i == 0
    end
    return tau.pow!(2).mul!(0.25)
  end
  
  def thermal_timescale(profile_id)
    data_filename = @data_dir + '/profile'+sprintf('%0.4d',profile_id)
    read_data(data_filename)
    t.do_box_labels(nil,'$P/g$ [\columnunit]','$\tau/\unitstyle{d}$')
    tau = eval_thermal_timescale
    ys = tau
    xs = @pressure.div(@gravity)
    t.xaxis_log_values = t.yaxis_log_values = true
    t.show_plot(set_boundaries(xs,ys,@margin)) do
      t.show_polyline(xs,ys)
    end
  end

  def eval_neutrino_timescale
    return @cp.mul(@temperature).div(@eps_nu)
  end
  
  def neutrino_timescale(profile_id)
    data_filename = @data_dir + '/profile'+sprintf('%0.4d',profile_id)
    read_data(data_filename)
    
    t.do_box_labels(nil,'$P/g$ [\columnunit]','$\tau_\nu/\unitstyle{d}$')
    
    xs = @pressure.div(@gravity).log10
    ys = eval_neutrino_timescale(profile_id)

    t.xaxis_log_values = t.yaxis_log_values = true
    t.show_plot(set_boundaries(xs,ys,@margin)) do
      t.show_polyline(xs,ys)      
    end
  end
  
  def show_both_timescales_with_legend(profile_id)
    t.show_plot_with_legend('plot_right_margin'=>0,'legend_left_margin'=>0.4,
      'legend_top_margin'=>0.6) { show_both_timescales(profile_id) }
  end
  
  def show_both_timescales(profile_id)
    data_filename = @data_dir + '/profile'+sprintf('%0.4d',profile_id)
    read_data(data_filename)
    t.do_box_labels(nil,'$P/g$ [\columnunit]','$\tau_\nu, \tau_{\mathrm{td}}$ [\unitstyle{d}]')
    xs = @pressure.div(@gravity).log10
    t.xaxis_log_values = t.yaxis_log_values = true
    ys_arry = [ eval_neutrino_timescale.div(86400.0).log10, eval_thermal_timescale.div(86400.0).log10 ]
    fmla = [
      '$\tau_\nu = CT/L_\nu$',
      '$\tau_{\mathrm{td}} = \frac{1}{4}\left[\int_r (\rho C/K)^{1/2}\dif r\right]^2$'
    ]
    t.show_plot(set_boundaries(xs,ys_arry[1],@margin)) do
        t.show_polyline(xs,ys_arry[0],Black,fmla[0],LINE_TYPE_SOLID)
        t.show_polyline(xs,ys_arry[1],Black,fmla[1],LINE_TYPE_DOTS)
    end
  end

  def density
  read_data
  t.do_box_labels(nil,'$P/g$ [\columnunit]','$\rho/(\grampercc)$')
  ys = @density.log10
  xs = @pressure.div(@grav).log10
  t.xaxis_log_values=t.yaxis_log_values=true
  t.show_plot(set_boundaries(xs,ys,@margin)) do
  t.show_polyline(xs,ys)
  end
  end

  def specific_heat
  read_data
  t.do_box_labels(nil,'$P/g$ [\columnunit]','$C_P\langle A\rangle/(\kB\NA)$')
  ys = @cp.mul(@abar).div(1.3806488e-16*6.02214129e23)
  xs = @pressure.div(@grav).log10
  t.xaxis_log_values=true
  t.show_plot(set_boundaries(xs,ys,@margin)) do
  t.show_polyline(xs,ys)
  end
  end

  def conductivity
  read_data
  t.do_box_labels(nil,'$P/g$ [\columnunit]','conductivity $K/(\ergs\nsp\second^{-1}\cm^{-2}\nsp\K^{-1})$')
  ys = @kcond.log10
  xs = @column.log10
  t.xaxis_log_values=t.yaxis_log_values=true
  t.show_plot(set_boundaries(xs,ys,@margin)) do
  t.show_polyline(xs,ys)
  end
  end

  def comp_with_zones
  read_data
  t.do_box_labels(nil,'$P/g$ ($10^{14}\nsp\columnunit$)','$\langle Z\rangle$, $\langle A\rangle$')
  xs = @pressure.div(@grav).div(1.0e14)
  zs = @zbar
  as = @abar
  zns = @column.div(1.0e14)
  t.xaxis_log_values = false
  t.show_plot(set_boundaries([0.5,3.0],as,@margin,ymin=20,ymax=56)) do
  t.context do
  t.stroke_color=LightSlateGrey
  t.stroke_width=0.5
  zns.each do |z|
  t.stroke_line(z,t.bounds_bottom,z,t.bounds_top)
  end
  end
  t.show_polyline(xs,zs,Red,nil,LINE_TYPE_SOLID)
  t.show_polyline(xs,as,Blue,nil,LINE_TYPE_SOLID)
  end
  end

  def enu_with_zones
  read_data
  t.do_box_labels(nil,'$P/g$ [\columnunit]','$\varepsilon_\nu/(\ergspersecond\nsp\gram^{-1})$')
  xs = @pressure.div(@grav).log10
  ys = @eps_nu.log10
  zns = @column.log10
  t.xaxis_log_values = true
  t.yaxis_log_values = true
  t.show_plot(set_boundaries([13.0,15.0],ys,@margin)) do
  t.context do
  t.stroke_color=LightSlateGrey
  t.stroke_width=0.5
  zns.each2(@zone) do |z,iz|
  t.stroke_line(z,t.bounds_bottom,z,t.bounds_top) unless iz % 2 != 0
  end
  end
  t.show_polyline(xs,ys)
  end
  end

  end

ProfilePlots.new
