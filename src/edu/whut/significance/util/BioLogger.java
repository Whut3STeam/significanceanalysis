package edu.whut.significance.util;

import java.io.File;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.logging.*;

/**
 * Created by Justin on 2014/12/21.
 */
public class BioLogger {
    private Logger m_log;
    private final boolean enableLoggerFile = false;

    public BioLogger(String path, String filename) {
        m_log = Logger.getLogger("significanceAnalysis");
//        m_log = Logger.getGlobal();
        m_log.setLevel(Level.ALL);

        ConsoleHandler consoleHandler = new ConsoleHandler();
        consoleHandler.setLevel(Level.ALL);
        consoleHandler.setFormatter(new BioLogHander());
        m_log.addHandler(consoleHandler);

        if (enableLoggerFile){
            try {
                String logfilename;
                if (!path.endsWith(File.separator)) {
                    logfilename = String.format("%s%s%s", path, File.separator, filename);
                } else {
                    logfilename = String.format("%s%s", path, filename);
                }
                //String logfilename = "/bai/ruanjun/Data/result.log";
                System.out.println(logfilename);
                FileHandler fileHandler = new FileHandler(logfilename, 500000, 10, true);

                fileHandler.setLevel(Level.ALL);

                fileHandler.setFormatter(new BioLogHander());
                m_log.addHandler(fileHandler);
            } catch (IOException e) {
                e.printStackTrace();
                System.err.println("\n\n\n" + e.getClass().getName() + ": " + e.getMessage());
                System.exit(0);
            }
        }

        m_log.setUseParentHandlers(false);
    }

    private class BioLogHander extends Formatter {
        private final DateFormat format = new SimpleDateFormat("hh:mm:ss:SSS");

        @Override
        public String format(LogRecord record) {
            return String.format("[%8d][%s] %s-> %s\n",
                    record.getSequenceNumber(), format.format(new Date()), record.getLevel(), record.getMessage());

        }
    }
}
