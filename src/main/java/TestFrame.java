/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 * @author Дмитрий
 */

import org.apache.commons.math3.linear.ArrayRealVector;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.*;
import java.util.List;

import static java.lang.Math.*;

public class TestFrame extends JFrame
{

    private Solver solver;

    public TestFrame(Solver solver)
    {
        super("Выбор графика");
        createGUI();
        this.solver = solver;

//        final Plotter t = new Plotter("test");
//        int size = 1000;
//        double[] arg = new double[size], foo = new double[size];
//        for(int i = 0; i<size; i++){
//            arg[i] = i/100.0;
//            foo[i] = solver.foo.value(arg[i]);
//        }
//        t.addRealVector(new ArrayRealVector(arg), new ArrayRealVector(foo), "function");
//        t.plot();

    }

    private double C = 2, K = 0.01, R = 2, T0 = 0, T1 = 1, alpha = 0.01;

    public static void main(String[] args)
    {
        javax.swing.SwingUtilities.invokeLater(new Runnable()
        {
            public void run()
            {
                JFrame.setDefaultLookAndFeelDecorated(false);
                TestFrame frame = new TestFrame(new Solver(
                        2,
                        0.01,
                        2,
                        0,
                        1,
                        0.01));
                frame.pack();
                frame.setLocationRelativeTo(null);
                frame.setVisible(true);
//                final Plotter plotter = new Plotter("tset", 640, 480);
//                frame.testC(10000, plotter);
//                plotter.plot();
            }
        });

        Solver solver = new Solver(
                2,
                0.01,
                2,
                0,
                1,
                0.01);
//        solver.testEps(Math.pow(10, -6), Math.pow(10,-9),Math.pow(10, -9));
        solver.testEps(Math.pow(10, -8), 0, pow(10, -4));
        solver.testEps(Math.pow(10, -8), 1.5, pow(10,-4));

//        Solver solver = new Solver(
//                        2,
//                        0.01,
//                        2,
//                        0,
//                        1,
//                        0.01);
//
//        solver.u(1,1,solver.getEigens(3000));

    }

    private int pointNum = 2000;
    private int time = 1000, radius = 2;
    private double[] arg, foo;
    double[] eigens;

    private List<Plotter> plotters = new LinkedList<>();

    public void plotAll()
    {
        for (Plotter p : plotters)
        {
            p.plot();
        }
    }

    public void testC(int num, Plotter plotter){
        if (eigens == null)
        {
            eigens = solver.getEigens(num);
        }
        arg = new double[num];
        foo = new double[num];
        double[] foo2 = new double[num];
        double[] foo3 = new double[num];

        double sum = 0;
        double res = solver.u(0.00000001, 0, eigens);
        for (int i = 0; i < eigens.length; i++)
        {
            arg[i] = i;
            foo[i] = Math.abs(solver.C(eigens[i]));// * solver.R(0.000001, eigens[i]));
//            foo[i] = (solver.chi(eigens[i]));// * solver.R(0.000001, eigens[i]));
//            foo[i] = Math.pow(Math.sin(eigens[i]*R/Math.sqrt(K)),2)/(2/(R*R) - 2*alpha/R);// * solver.R(0.000001, eigens[i]));
//            foo[i] = (Math.sin(eigens[i])/Math.pow(eigens[i],2));// * solver.R(0.000001, eigens[i]));
//            System.out.println(foo[i]);
//            foo[i] = Math.abs(solver.u(0.0001, 0, eigens));
//            foo2[i] = Math.abs(solver.R(0.1, eigens[i]));
            foo2[i] = Math.abs(solver.R(pow(10, -9), eigens[i]));
//            foo2[i] = 8*(Math.sqrt(1 + Math.pow(Math.PI*(i+1),2)))/(2*Math.pow(Math.PI*i,2) - i);
            foo3[i] = (solver.T(0, eigens[i]));
//            System.out.println(foo[i]*foo2[i]*foo3[i]);
//            System.out.println(solver.chi(eigens[i]));
            sum += Math.pow(-1,i)*foo[i]*foo2[i]*foo3[i];
//            System.out.println(sum -1 - res);
        }
        plotter.addRealVector(new ArrayRealVector(arg), new ArrayRealVector(foo), "C=" + num);
//        plotter.addRealVector(new ArrayRealVector(arg), new ArrayRealVector(foo2), "R=" + num);
//        plotter.addRealVector(new ArrayRealVector(arg), new ArrayRealVector(foo3), "T=" + num);
    }

    public void constantRadius(double rad, Plotter plotter)
    {
        if (eigens == null)
        {
            eigens = solver.getEigens(pointNum);
        }
        arg = new double[time];
        foo = new double[time];

        for (int i = 0; i < time; i++)
        {
            arg[i] = i;
            foo[i] = solver.u(rad, i, eigens);
        }
        System.out.println("finished calc for rad=" + rad);
        plotter.addRealVector(new ArrayRealVector(arg), new ArrayRealVector(foo), "Радиус=" + rad);
    }

    public static final double TIME_STEP = pow(10, -2);

    public void constantTime(double t, Plotter plotter)
    {
        System.out.println("constantTime start");
        if (eigens == null)
        {
            eigens = solver.getEigens(pointNum);
        }
        System.out.println("eigens ready");
        int size = (int) (R / TIME_STEP);
        arg = new double[size];
        foo = new double[size];

        double x = TIME_STEP;
        for (int i = 0; i < size; i++)
        {
            arg[i] = x;
            foo[i] = solver.u(x, t, eigens);
            x += TIME_STEP;
        }
        System.out.println("finished calc for tune =" + t);
        plotter.addRealVector(new ArrayRealVector(arg), new ArrayRealVector(foo), "Время=" + t);
    }

    private JButton addButton(String label, JPanel panel)
    {
        final JButton btn = new JButton(label);
        panel.add(btn);
        return btn;
    }

    private JTextArea addText(String text, JPanel panel)
    {
        final JTextArea area = new JTextArea(text);
        area.setSize(100, 30);
        panel.add(area);
        return area;
    }

    private JTextField addTextField(String text, JPanel panel)
    {
        final JTextField field = new JTextField(text);
        field.setSize(100, 30);
        panel.add(field);
        return field;
    }

    JTextField inpK, inpC, inpAlpha, inpT0, inpT1, inpRadius, inpR, inpT;

    public void createGUI()
    {
//        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        JPanel btnPanel = new JPanel();
        btnPanel.setLayout(new BoxLayout(btnPanel, BoxLayout.Y_AXIS));

        JButton btnBuildTime = addButton("Построить график для заданного времени", btnPanel);
        btnBuildTime.addMouseListener(new MouseAdapter()
        {
            @Override
            public void mouseClicked(MouseEvent e)
            {
                double t = 1;
                try
                {
                    t = Double.parseDouble(inpT.getText());
                } catch (Exception ex)
                {
                    JOptionPane.showMessageDialog(TestFrame.this, "Ошибка при вводе данных");
                }
                final Plotter plotter = new Plotter("Время=" + t, 640, 480);
                constantTime(t, plotter);
                plotter.plot();
            }
        });

        JButton btnBuildRadius = addButton("Построить график для заданного радиуса", btnPanel);
        btnBuildRadius.addMouseListener(new MouseAdapter()
        {
            @Override
            public void mouseClicked(MouseEvent e)
            {
                double r = 1;
                try
                {
                    r = Double.parseDouble(inpR.getText());
                } catch (Exception ex)
                {
                    JOptionPane.showMessageDialog(TestFrame.this, "Ошибка при вводе данных");
                }
                final Plotter plotter = new Plotter("Радиус=" + r, 640, 480);
                constantRadius(r, plotter);
                plotter.plot();
            }
        });

        JButton applyData = addButton("Применить данные", btnPanel);
        applyData.addMouseListener(new MouseAdapter()
        {
            @Override
            public void mouseClicked(MouseEvent e)
            {
                try
                {
                    solver = new Solver(
                            Double.parseDouble(inpC.getText()),
                            Double.parseDouble(inpK.getText()),
                            Double.parseDouble(inpRadius.getText()),
                            Double.parseDouble(inpT0.getText()),
                            Double.parseDouble(inpT1.getText()),
                            Double.parseDouble(inpAlpha.getText())
                    );

                } catch (Exception ex)
                {
                    JOptionPane.showMessageDialog(TestFrame.this, "Ошибка ввода");
                }
            }
        });

        getContentPane().add(btnPanel);

        JPanel textPanel = new JPanel();
        textPanel.setLayout(new BoxLayout(textPanel, BoxLayout.Y_AXIS));

        JTextArea msgK = addText("Параметр К", textPanel),
                msgC = addText("Параметр С", textPanel),
                msgAlpha = addText("Параметр альфа", textPanel),
                msgT0 = addText("Начальная температура", textPanel),
                msgT1 = addText("Конечная температура", textPanel),
                msgRadius = addText("Радиус", textPanel),
                msgR = addText("Радиус для расчета", textPanel),
                msgT = addText("Время для расчета", textPanel);
        getContentPane().add(textPanel);
        textPanel.setPreferredSize(new Dimension(100, 200));

        JPanel inpPanel = new JPanel();
        inpPanel.setLayout(new BoxLayout(inpPanel, BoxLayout.Y_AXIS));

        inpK = addTextField(String.valueOf(K), inpPanel);
        inpC = addTextField(String.valueOf(C), inpPanel);
        inpAlpha = addTextField(String.valueOf(alpha), inpPanel);
        inpT0 = addTextField(String.valueOf(T0), inpPanel);
        inpT1 = addTextField(String.valueOf(T1), inpPanel);
        inpRadius = addTextField(String.valueOf(radius), inpPanel);
        inpR = addTextField("", inpPanel);
        inpT = addTextField("", inpPanel);
        inpPanel.setPreferredSize(new Dimension(100, 200));

        getContentPane().add(inpPanel);

        setLayout(new FlowLayout());

        setPreferredSize(new Dimension(400, 800));
    }
}