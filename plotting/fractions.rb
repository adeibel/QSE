# plots.rb

require 'Tioga/FigureMaker'
require './plot_styles.rb'

class MyPlots

    include Math
    include Tioga
    include FigureConstants
    include MyPlotStyles
    
    def t
        @figure_maker
    end

    def initialize
        @figure_maker = FigureMaker.default

        t.save_dir = 'plots_out'
        t.def_eval_function { |str| eval(str) }
        
        @data_filename = "y_output.data"
       # @data_filename2 = "ashes_outer_crust.data"
		@data_filename2 = "outer_crust.txt"
        @data_filename3 = "ashes_acc.txt"
        @data_filename4 = "SLy4_0_14_freq.o2.txt"
        @opacity_data = nil
        
        @image_right_margin = 0.07
        @margin = 0.0
        @header = Dvector.new
        @positions = Dvector.new
        @blues = Dvector.new
        @reds = Dvector.new
        @greens = Dvector.new
        @data_array = [@positions, @blues, @reds, @greens]
        @have_data = false
        @big_blues = Dvector.new
        @big_blues_scale = 11.0
	    @pressure = Dvector.new
	    @nb = Dvector.new
	    @Ye = Dvector.new
	    @Yn = Dvector.new
	    @Z = Dvector.new
	    @A = Dvector.new
	    @mue = Dvector.new
	    @mun = Dvector.new
	    @crust_data = [@pressure, @nb, @Ye, @Yn, @Z, @A, @mue, @mun]
	    
	    @pressure2 = Dvector.new
	    @Z_outer = Dvector.new
	    @A_outer = Dvector.new
	    @outer_crust = [@pressure2, @Z_outer, @A_outer]
	    
	    @pressure3 = Dvector.new
	    @n_b2 = Dvector.new
	    @Ye2 = Dvector.new
	    @Yn2 = Dvector.new
	    @fr = Dvector.new
	    @Z_equil = Dvector.new
	    @A_equil = Dvector.new
	    @equil_outer = [@pressure3, @n_b2, @Ye2, @Yn2, @fr, @Z_equil, @A_equil]
 
	    @pressure4 = Dvector.new
	    @nb3 = Dvector.new
	    @N2 = Dvector.new
	    @Z2 = Dvector.new
	    @equil = [@pressure4, @nb3, @N2, @Z2]
              
        t.def_figure("fractions") { plot1 }

        t.model_number = -1
        
        t.def_enter_page_function { enter_page }

        # Uncomment the followiing line if you want tioga to leave
        # its temporary files behind.
        t.autocleanup = false
            
    end
    
    def enter_page
        set_default_plot_style
        t.default_page_width = 72*3 # in big-points (1/72 inch)
        t.default_page_height = t.default_page_width
        t.default_enter_page_function
    end

    def french_decimal_separator
      t.tex_preamble += <<'EOD'
\def\frenchsep#1{\frenchsepma#1\endoffrenchsep} 
\def\fseat#1{\frenchsepma}
\def\frenchsepma{\futurelet\next\frenchsepmw}
\def\frenchsepmw{\ifx\next\endoffrenchsep\let\endoffrenchsep=\relax%
\else\if\next.\ifmmode\mathord,\else,\fi%
\else\next\fi\expandafter
\fseat\fi}
EOD

      t.yaxis_numeric_label_tex = '$\frenchsep{#1}$'

      #blues
    plot1
    end
    
    def read_data
        Dvector.read(@data_filename, @crust_data, 2)
        Dvector.read(@data_filename2, @outer_crust,1)
        Dvector.read(@data_filename3, @equil_outer, 3)
        Dvector.read(@data_filename4, @equil,3)
        @have_data = true
        t.need_to_reload_data = false
    end

#     def show_model_number(pos = 1, shift = 4)
#         if !(t.in_subplot) && t.model_number > 0
#             t.show_text('text' => t.model_number.to_s,
#                 'side' => TOP, 'pos' => pos,
#                 'shift' => shift, 'scale' => 0.8,
#                 'justification' => RIGHT_JUSTIFIED)
#         end
#     end
    
  def plot_boundaries(xs, ys, ymin=nil, ymax=nil)
    xmin = -4.0 #xs.min
    xmax = xs.max
    ymin = 0.0 #ys.min if ymin == nil
    ymax = 1.0 #ys.max if ymax == nil
    width = (xmax == xmin)? 1 : xmax - xmin
    height = (ymax == ymin)? 1 : ymax - ymin
    left_boundary = xmin #-1.5
    right_boundary = xmax #+0.5
    top_boundary = ymax #+0.5
    bottom_boundary = ymin #-1.5
    return [left_boundary,right_boundary,top_boundary,bottom_boundary]
  end

    
    def plot1
        read_data
        t.do_box_labels('', 'pressure [MeV fm$^{-3}$]', 'fraction')
        xs_in = @pressure
        xs = (xs_in).log10
        ys = @Z
        ys2 = @A
        ys_n = (@A).minus(@Z)
        ye1 = @Ye
        yn1 = @Yn
        
        xs_outer = @pressure2
        xs2 = (xs_outer).log10
        ys3 = @Z_outer
        ys4 = @A_outer
        ye2 = @Ye2
        yn2 = @Yn2
        
        xs_outer2 = @pressure3
        xs3 = (xs_outer2).log10
        ys5 = @Z_equil
        ys6 = @A_equil
        
        xs_equil = @pressure4
        xs4 = (xs_equil).log10
        ys7 = @Z2
        ys8 = @Z2+@N2
        
        t.xaxis_major_tick_length = 0.25
        t.xaxis_minor_tick_length = 0.15
        t.yaxis_major_tick_length = 0.25
        t.yaxis_minor_tick_length = 0.15
        
        t.xaxis_log_values = true
		t.line_width = 0.5

        t.show_plot(plot_boundaries(xs,ys,@margin)){
		t.show_polyline(xs,ye1,Black)
		t.show_polyline(xs,yn1,Blue)
        }
        
#         t.show_label('text' => '$<$Z$>$', 'x' => -5.5 , 'y' => 45, 
#             'justification' => CENTERED, 'color' => Black, 'scale' => 1.0) 
#         t.show_label('text' => '$<$A$>$', 'x' => -5.5 , 'y' => 95, 
#             'justification' => CENTERED, 'color' => Black, 'scale' => 1.0)
    end
    

end



MyPlots.new
